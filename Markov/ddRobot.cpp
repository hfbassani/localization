/**
 * Robô de direção diferencial
 * Disciplina de Robótica CIn/UFPE
 * 
 * @autor Prof. Hansenclever Bassani
 * 
 * Este código é proporcionado para facilitar os passos iniciais da programação.
 * Porém, não há garantia de seu correto funcionamento.
 * 
 * Testado em: Ubuntu 14.04 + Netbeans
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <unistd.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

#define V_REP_IP_ADDRESS "10.0.2.2"//"127.0.0.1"
#define V_REP_PORT 19997//1999;

extern "C" {
#include "extApi.h"
    /*	#include "extApiCustom.h" if you wanna use custom remote API functions! */
}

using namespace std;

simxInt ddRobotHandle;
simxInt leftMotorHandle;
simxInt rightMotorHandle;
simxInt graphOdometryHandle;

simxInt sensorFrontHandle;
simxInt sensorLeftHandle;
simxInt sensorRightHandle;

float sensorFrontPos[3] = {0.03, 0, 0};
float sensorLeftPos[3]  = {0.03, 0, M_PI/2};
float sensorRightPos[3] = {0.03, 0, -M_PI/2};
        
#define BEL_NXY 40
#define BEL_NTHETA (360/8)
#define SIGMA 0.2
#define GSCALE  0.089422044
#define GMSCALE 0.089422044

//Robot parameters:
#define r 0.02
#define l 0.1

#define Kr 0.1
#define Kl 0.1

#define MINPROB 0.001

float mapLines[20][4];
int nLines=0;

cv::Scalar lineColor(255, 0, 0);
cv::Scalar pointColor(0,0,255);
        
void getPosition(int clientID, simxFloat pos[]) { //[x,y,theta]

    simxInt ret = simxGetObjectPosition(clientID, ddRobotHandle, -1, pos, simx_opmode_oneshot_wait);
    if (ret > 0) {
        printf("Error reading robot position\n");
        return;
    }

    simxFloat orientation[3];
    ret = simxGetObjectOrientation(clientID, ddRobotHandle, -1, orientation, simx_opmode_oneshot_wait);
    if (ret > 0) {
        printf("Error reading robot orientation\n");
        return;
    }

    simxFloat theta = orientation[2];
    pos[2] = theta;
}

simxInt getSimTimeMs(int clientID) { //In Miliseconds
    return simxGetLastCmdTime(clientID);
}

float to_positive_angle(float angle) {

    angle = fmod(angle, 2 * M_PI);
    while (angle < 0) {
        angle = angle + 2 * M_PI;
    }
    return angle;
}

float smallestAngleDiff(float target, float source) {
    float a;
    a = to_positive_angle(target) - to_positive_angle(source);

    if (a > M_PI) {
        a = a - 2 * M_PI;
    } else if (a < -M_PI) {
        a = a + 2 * M_PI;
    }
    return a;
}

void readOdometers(int clientID, simxFloat &dwL, simxFloat &dwR) {
    //old joint angle position
    static simxFloat lwprev=0; 
    static simxFloat rwprev=0;
    
    //current joint angle position
    simxFloat lwcur=0;
    simxFloat rwcur=0;

    simxGetJointPosition(clientID, leftMotorHandle, &lwcur, simx_opmode_oneshot);
    simxGetJointPosition(clientID, rightMotorHandle, &rwcur, simx_opmode_oneshot);

    dwL = smallestAngleDiff(lwcur, lwprev);
    dwR = smallestAngleDiff(rwcur, rwprev);
    lwprev = lwcur;
    rwprev = rwcur;
}

simxFloat readSonar(int clientID, simxInt sensorHandle) {    
    simxFloat detectedPoint[3];
    simxUChar detectionState=1;
    simxInt ret = simxReadProximitySensor(clientID, sensorHandle, &detectionState, detectedPoint, NULL, NULL, simx_opmode_buffer);
    
    //printf("ret: %d ds: %d\n", ret, detectionState);
    if (ret<=0 && detectionState==1)
        return detectedPoint[2];
 
    return -1;
}

void setTargetSpeed(int clientID, simxFloat phiL, simxFloat phiR) {
    simxSetJointTargetVelocity(clientID, leftMotorHandle, phiL, simx_opmode_oneshot);
    simxSetJointTargetVelocity(clientID, rightMotorHandle, phiR, simx_opmode_oneshot);   
}

inline double to_deg(double radians) {
    return radians * (180.0 / M_PI);
}

void sensorToRobot(float dist, float* sensorPos, float *out) {
    out[0] = sensorPos[0] + dist*cos(sensorPos[2]);
    out[1] = sensorPos[1] + dist*sin(sensorPos[2]);
}

void robotToWorld(float* vetin, float* robotPos, float *out) {
    float sintheta = sin(robotPos[2]);
    float costheta = cos(robotPos[2]);
    
    out[0] = robotPos[0] + (vetin[0]*costheta - vetin[1]*sintheta);
    out[1] = robotPos[1] + (vetin[0]*sintheta + vetin[1]*costheta);
    out[2] = robotPos[2] + vetin[2];
}

void worldToMap(float *vetin, float *out) {
    
    out[0] = 50+(vetin[0]+2)*400/4;
    out[1] = 50+(2 - vetin[1])*400/4;
}

float distance2(float *p1, float* p2) {
    return (p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]);
}

float distance(float *p1, float* p2) {
    return sqrt((p1[0]-p2[0])*(p1[0]-p2[0]) + (p1[1]-p2[1])*(p1[1]-p2[1]));
}

int nearesPointInLine(float *p, float* pA, float* pB, float *npl) {
  float AP[2] = {p[0] - pA[0], p[1] - pA[1]}; // Storing vector A->P
  float AB[2] = {pB[0] - pA[0], pB[1] - pA[1]}; // Storing vector A->B

  float ab2 = AB[0]*AB[0] + AB[1]*AB[1];  //Basically finding the squared magnitude of AB

  float AP_dot_AB = AP[0]*AB[0] + AP[1]*AB[1]; // The dot product of AP and AB                                     

  float t = AP_dot_AB / ab2; // The normalized "distance" from a to
  
  npl[0] = pA[0] + AB[0]*t;
  npl[1] = pA[1] + AB[1]*t;
  
//  printf("p: [%.2f, %.2f]\n", p[0], p[1]);
//  printf("pA: [%.2f, %.2f]\n", pA[0], pA[1]);
//  printf("pB: [%.2f, %.2f]\n", pB[0], pB[1]);
//  printf("AP: [%.2f, %.2f]\n", AP[0], AP[1]);
//  printf("AB: [%.2f, %.2f]\n", AB[0], AB[1]);
//  printf("npl: [%.2f, %.2f]\n", npl[0], npl[1]);

  if (t<0 || t>sqrt(ab2))
      return 0; //point is not inside AB

  return 1; //point is inside AB
}

float nearesPointInMap(float *p, float *npl) {
    float pl[2], d, best_d=99999;
    
    //printf("\n");
    for (int i=0; i<nLines; i++) {
        nearesPointInLine(p, &mapLines[i][0], &mapLines[i][2], pl);
        {
            d = distance2(pl, p);
            if (d<best_d) {
                best_d = d;
                npl[0] = pl[0];
                npl[1] = pl[1];
            }
//            printf("p: [%.2f, %.2f] ", p[0], p[1]);
//            printf("pl: [%.2f, %.2f] ", pl[0], pl[1]);
//            printf("d:%.4f bd: %.4f ", d, best_d);
//            printf("npl: [%.2f, %.2f]\n", npl[0], npl[1]);
        }
    }
    
    return best_d;
}

void robotPosFromWall(float d, float *npl, float *pRobot, float *newRobot) {
    float NPLR[2] = {pRobot[0] - npl[0], pRobot[1] - npl[1]};
    
    float norm = sqrt(NPLR[0]*NPLR[0] + NPLR[1]*NPLR[1]);
    newRobot[0] = npl[0] + d*NPLR[0]/norm;
    newRobot[1] = npl[1] + d*NPLR[1]/norm;
    
    //printf("NPLR: [%.2f, %.2f]\n", NPLR[0], NPLR[1]);
    //printf("newPos: [%.2f, %.2f]\n", newRobot[0], newRobot[1]);
}

void plotData(cv::Mat &image, float* pointW) {
    float pointM[2];
    worldToMap(pointW, pointM);
    
    if (pointM[0]>=0 && pointM[0]<image.rows && pointM[1]>=0 && pointM[1]<image.cols) {
        cv::circle(image, cv::Point(pointM[0],pointM[1]), 2, pointColor, 1);
    } else
        printf("Invalid pointM: [%.2f, %.2f]\n", pointM[0], pointM[1]);
}

void plotLine(cv::Mat &image, float* pointWA, float* pointWB) {
    float pointMA[2], pointMB[2];
    worldToMap(pointWA, pointMA);
    worldToMap(pointWB, pointMB);
    
//    printf("LineA: [%.2f, %.2f]\n", pointMA[0], pointMA[1]);
//    printf("LineB: [%.2f, %.2f]\n", pointMB[0], pointMB[1]);
    
    cv::Point pa(pointMA[0],pointMA[1]);
    cv::Point pb(pointMB[0],pointMB[1]);
    
    cv::line(image, pa, pb, lineColor, 2);
}

int loadMap(const char *filename) {

    FILE *fp = fopen(filename, "r");
    if (fp !=NULL) {
        nLines = 0;
        while (1) {
            fscanf(fp, "%f%f%f%f", &mapLines[nLines][0], &mapLines[nLines][1], &mapLines[nLines][2], &mapLines[nLines][3]);
            if (feof(fp)) break;
            printf("line loaded: (%.2f, %.2f) (%.2f, %.2f)\n", mapLines[nLines][0], mapLines[nLines][1], mapLines[nLines][2], mapLines[nLines][3]);
            nLines++;
        }
        fclose(fp);
        return 1;
    } else
        printf("Could not open file: %s\n", filename);

    return 0;
}

int drawMap(cv::Mat &image) {

    for (int i=0; i<nLines; i++) {
        plotLine(image, &mapLines[i][0], &mapLines[i][2]);
    }
}

void robotToSensorPoint(float *pRobot, float *pSensor, float dist, float* point){
    float pSR[2];
    sensorToRobot(dist, pSensor, pSR);
    robotToWorld(pSR, pRobot, point);
}

float pGaussian(float dist) {
    return GSCALE*exp(-dist*dist/(2*SIGMA*SIGMA))/(SIGMA*sqrt(2*M_PI));    
}

float pMultGaussian(cv::Mat &m, cv::Mat &x, cv::Mat E) {
    cv::Mat x_m = x-m;
    cv::Mat Ei = E.inv();
    
//    cout << "m: (" << m.rows << "," << m.cols << ") x: (" << x.rows << "," << x.cols << ") E: (" << E.rows << "," << E.cols << ")" << endl;
//    cout << "x-m: (" << x_m.rows << "," << x_m.cols << ") Ei: (" << Ei.rows << "," << Ei.cols << ") " << endl;
    
    cv::Mat D = (x_m*Ei)*x_m.t();//*x_m;
    float e = exp(-D.at<float>(0,0)/2);
    float det = cv::determinant(E);
    
    return (GMSCALE*e/(pow(2*M_PI,x.cols/2.0)*sqrt(det)));
}

void perceptionUpdate(cv::Mat &bel, float distF, float distL, float distR){
    float robotPos[3]; //[x,y,theta]
    int x,y,t;
    float sensorPoint[3], mapPoint[3];
    float dF, dL, dR, p, sum=0;
    
    for(robotPos[0]=-2, x=0; x<BEL_NXY; robotPos[0]+=4.0/BEL_NXY, x++)
    for(robotPos[1]=-2, y=0; y<BEL_NXY; robotPos[1]+=4.0/BEL_NXY, y++)
    for(robotPos[2]=-M_PI, t=0; t<BEL_NTHETA; robotPos[2]+=2*M_PI/BEL_NTHETA, t++) {
        
        robotToSensorPoint(robotPos, sensorFrontPos, distF, sensorPoint);
        dF = nearesPointInMap(sensorPoint, mapPoint);
        
        robotToSensorPoint(robotPos, sensorLeftPos, distL, sensorPoint);
        dL = nearesPointInMap(sensorPoint, mapPoint);
        
        robotToSensorPoint(robotPos, sensorRightPos, distR, sensorPoint);
        dR = nearesPointInMap(sensorPoint, mapPoint);
        
        //printf("d: %f, p(d): %f\n", dF, pGaussian(dF));
        p = pGaussian(dF)*pGaussian(dL)*pGaussian(dR)*bel.at<float>(x,y,t);
        bel.at<float>(x,y,t) = p;
        sum+=p;
    }        
    
    bel = bel/sum;
}


void fdeltaRL(float theta, float ds, float dtheta, cv::Mat &FDrl) {
    float costdt2 = cos(theta+dtheta/2);
    float sintdt2 = sin(theta+dtheta/2);
    float b = 2*l;
    
    FDrl.at<float>(0,0) = costdt2/2 - (ds/(2*b))*sintdt2;//1
    FDrl.at<float>(0,1) = costdt2/2 + (ds/(2*b))*sintdt2;//2
    
    FDrl.at<float>(1,0) = sintdt2/2 + (ds/(2*b))*costdt2;//3
    FDrl.at<float>(1,1) = sintdt2/2 - (ds/(2*b))*costdt2;//4

    FDrl.at<float>(2,0) = 1/b;//5
    FDrl.at<float>(2,1) = -1/b;//6
    

}

void odomError(cv::Mat &FDrl, cv::Mat &Ed, cv::Mat &Ep) {
    cv::Mat FE = FDrl*Ed;
    cv::Mat FDrlT = FDrl.t();
    Ep = FE*FDrlT;
}

void actionUpdate(cv::Mat &bel, float dsl, float dsr) {
    int x1,y1,t1; //x,y,theta indexes on bel map
    int x0,y0,t0; //x,y,theta indexes on tempbel map
    
    float dtheta = (dsr - dsl)/(2*l);
    float ds = (dsr + dsl)/2;
    cv::Mat X1(1,3, CV_32FC1), //current position [x,y,theta]
            X0(1,3, CV_32FC1), //previous position [x,y,theta]
            M(1,3, CV_32FC1),  //new expected position [x,y,theta]
            Ed(2,2, CV_32FC1),  //covar(dsr, dsl)
            FDrl(3,2, CV_32FC1),//Jacobian
            Ep(3,3, CV_32FC1);  //motion error Sigmap
    
    const int belSizes[3]={BEL_NXY, BEL_NXY, BEL_NTHETA};
    cv::Mat sumbel(3, belSizes, CV_32FC1, 0.0); //somatory belief mat init with zeroes
        
    Ed.at<float>(0,0) = Kr*dsr;
    Ed.at<float>(1,1) = Kl*dsl;
    Ed.at<float>(0,1) = 0;
    Ed.at<float>(1,0) = 0;
    
//    for(X1.at<float>(0,2)=-M_PI, t=0; t<BEL_NTHETA; X1.at<float>(0,2)+=2*M_PI/BEL_NTHETA, t++) {//for each new theta
//        
//        float costdt2 = cos(X1.at<float>(0,2)+dtheta/2);
//        float sintdt2 = sin(X1.at<float>(0,2)+dtheta/2);
//        fdeltaRL(X1.at<float>(0,2), ds, dtheta, FDrl);
//        odomError(FDrl, Ed, Ep);
//        
//        for(X1.at<float>(0,0)=-2, x=0; x<BEL_NXY; X1.at<float>(0,0)+=4.0/BEL_NXY, x++)                         
//        for(X1.at<float>(0,1)=-2, y=0; y<BEL_NXY; X1.at<float>(0,1)+=4.0/BEL_NXY, y++) { //for each new (x,y)
//                
//            for(X0.at<float>(0,2)=-M_PI, t0=0; t0<BEL_NTHETA; X0.at<float>(0,2)+=2*M_PI/BEL_NTHETA, t0++) {//for each old theta
//                M.at<float>(0,2) = X0.at<float>(0,2)+dtheta;
//
//                for(X0.at<float>(0,0)=-2, x0=0; x0<BEL_NXY; X0.at<float>(0,0)+=4.0/BEL_NXY, x0++) { //for each old x
//                    M.at<float>(0,0) = X0.at<float>(0,0)+(ds/2)*costdt2;
//
//                    for(X0.at<float>(0,1)=-2, y0=0; y0<BEL_NXY; X0.at<float>(0,1)+=4.0/BEL_NXY, y0++) {  //for each old y
//                        M.at<float>(0,1) = X0.at<float>(0,1)+(ds/2)*sintdt2;
//
//                        float px1_u1x0 = pMultGaussian(M, X1, Ep);
//                        //sumbel.at<float>(x0,y0,t0) += bel.at<float>(x0,y0,t0)*px1_u1x0;
//                    }
//                }
//            }
//        }
//    }
    
    float sum=0;
    for(X0.at<float>(0,2)=-M_PI, t0=0; t0<BEL_NTHETA; X0.at<float>(0,2)+=2*M_PI/BEL_NTHETA, t0++) {//for each old theta
        M.at<float>(0,2) = X0.at<float>(0,2)+dtheta;
        
        for(X0.at<float>(0,0)=-2, x0=0; x0<BEL_NXY; X0.at<float>(0,0)+=4.0/BEL_NXY, x0++) { //for each old x
            for(X0.at<float>(0,1)=-2, y0=0; y0<BEL_NXY; X0.at<float>(0,1)+=4.0/BEL_NXY, y0++) {  //for each old y

                float b = bel.at<float>(x0,y0,t0);
                if (b>MINPROB)
                for(X1.at<float>(0,2)=-M_PI, t1=0; t1<BEL_NTHETA; X1.at<float>(0,2)+=2*M_PI/BEL_NTHETA, t1++) {//for each new theta
                            
                    float costdt2 = cos(X1.at<float>(0,2)+dtheta/2);
                    float sintdt2 = sin(X1.at<float>(0,2)+dtheta/2);
                    fdeltaRL(X1.at<float>(0,2), ds, dtheta, FDrl);
                    odomError(FDrl, Ed, Ep);                    

                    for(X1.at<float>(0,0)=-2, x1=0; x1<BEL_NXY; X1.at<float>(0,0)+=4.0/BEL_NXY, x1++) { //for each new x
                        M.at<float>(0,0) = X0.at<float>(0,0)+(ds/2)*costdt2;

                        for(X1.at<float>(0,1)=-2, y1=0; y1<BEL_NXY; X1.at<float>(0,1)+=4.0/BEL_NXY, y1++) { //for each new y           
                            M.at<float>(0,1) = X0.at<float>(0,1)+(ds/2)*sintdt2;

                            float px1_u1x0 = pMultGaussian(M, X1, Ep);
                            cout << "px: " << px1_u1x0 << endl;
                            sumbel.at<float>(x1,y1,t1) += b*px1_u1x0;
                            sum += b*px1_u1x0;
                        }
                    }
                }
            }
        }
    }
    
    bel = sumbel/sum;
}

void drawBel(cv::Mat &image, cv::Mat bel) {
    float robotPos[3]; //[x,y,theta]
    float map[2];
    
    int x,y,t;
    float p, min=9999, max=0;
    
    for(x=0; x<BEL_NXY; x++)
    for(y=0; y<BEL_NXY; y++)
    for(t=0; t<BEL_NTHETA;  t++) {
        p=bel.at<float>(x, y, t);
        if (p>max) max = p;
        if (p<min) min = p;
    }
            
    for(robotPos[0]=-2, x=0; x<BEL_NXY; robotPos[0]+=4.0/BEL_NXY, x++)
    for(robotPos[1]=-2, y=0; y<BEL_NXY; robotPos[1]+=4.0/BEL_NXY, y++) {
        float p=0;
        for(t=0; t<BEL_NTHETA;  t++) {
            if (bel.at<float>(x, y, t)>p)
                p = bel.at<float>(x, y, t);
        }
        
        worldToMap(robotPos, map);
        int color = (p-min)/(max-min)*255;
        //printf("p: %f, min: %f, max: %f\n", p, min, max);
        //printf("p: (%f, %f)\n", map[0],map[1]);
        image.at<cv::Vec3b>(map[0],map[1])[0] = color;
        image.at<cv::Vec3b>(map[0],map[1])[1] = color;
        image.at<cv::Vec3b>(map[0],map[1])[2] = color;
    }
}

int main(int argc, char* argv[]) {

    std::string ipAddr = V_REP_IP_ADDRESS;
    int portNb = V_REP_PORT;

    if (argc > 1) {
        ipAddr = argv[1];
    }

    printf("Iniciando conexao com: %s...\n", ipAddr.c_str());

    int clientID = simxStart((simxChar*) (simxChar*) ipAddr.c_str(), portNb, true, true, 2000, 5);
    if (clientID != -1) {
        printf("Conexao efetuada\n");
        
        //Get handles for robot parts, actuators and sensores:
        simxGetObjectHandle(clientID, "RobotFrame#", &ddRobotHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "LeftMotor#", &leftMotorHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "RightMotor#", &rightMotorHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "GraphOdometry#", &graphOdometryHandle, simx_opmode_oneshot_wait);

        simxGetObjectHandle(clientID, "ProximitySensorF#", &sensorFrontHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "ProximitySensorL#", &sensorLeftHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "ProximitySensorR#", &sensorRightHandle, simx_opmode_oneshot_wait);
        
        printf("RobotFrame: %d\n", ddRobotHandle);
        printf("LeftMotor: %d\n", leftMotorHandle);
        printf("RightMotor: %d\n", rightMotorHandle);
        printf("GraphOdometry: %d\n", graphOdometryHandle);
        
        printf("ProximitySensorF: %d\n", sensorFrontHandle);
        printf("ProximitySensorL: %d\n", sensorLeftHandle);
        printf("ProximitySensorR: %d\n", sensorRightHandle);

        //start simulation
        int ret = simxStartSimulation(clientID, simx_opmode_oneshot_wait);
        
        if (ret==-1) {
            printf("Não foi possível iniciar a simulação.\n");
            return -1;
        }
        
        printf("Simulação iniciada.\n");

        simxFloat detectedPoint[3];
        simxUChar detectionState=1;
        simxReadProximitySensor(clientID, sensorFrontHandle, &detectionState, detectedPoint, NULL, NULL, simx_opmode_streaming);
        simxReadProximitySensor(clientID, sensorLeftHandle, &detectionState, detectedPoint, NULL, NULL, simx_opmode_streaming);
        simxReadProximitySensor(clientID, sensorRightHandle, &detectionState, detectedPoint, NULL, NULL, simx_opmode_streaming);

        cv::Mat image(500, 500, CV_8UC3, cv::Scalar(0,0,0));
        
        if (loadMap("map.txt"))
            drawMap(image);
        
        lineColor[0] = 0;
        
        const int belSizes[3]={BEL_NXY, BEL_NXY, BEL_NTHETA};
        cv::Mat bel(3, belSizes, CV_32FC1, (1.0/(BEL_NXY*BEL_NXY*BEL_NTHETA)));
        
        //While is connected:
        while (simxGetConnectionId(clientID) != -1)
        {  
            //Read current position:
            simxFloat pos[3]; //[x,y,theta] in [cm cm rad]
            getPosition(clientID, pos);

            //Read simulation time of the last command:
            simxInt time = getSimTimeMs(clientID); //Simulation time in ms or 0 if sim is not running
            //stop the loop if simulation is has been stopped:
            if (time == 0) break;             
            //printf("Posicao: [%.2f %.2f %.2fº], time: %dms\n", pos[0], pos[1], to_deg(pos[2]), time);
            
            //Read current wheels angle variation:
            simxFloat dPhiL, dPhiR; //rad
            readOdometers(clientID, dPhiL, dPhiR);
            //printf("dPhiL: %.2f dPhiR: %.2f\n", dPhiL, dPhiR);
            
            simxFloat distF=-1, distR=-1, distL=-1;
            distF = readSonar(clientID, sensorFrontHandle);
            distL = readSonar(clientID, sensorLeftHandle);
            distR = readSonar(clientID, sensorRightHandle);
            
            if (distF>=0 && distL>=0 && distR >=0) {
                printf("start perception update...");fflush(stdout);
                perceptionUpdate(bel, distF, distL, distR);
                printf("done\n");
            }
            
            printf("start action update...");fflush(stdout);
            actionUpdate(bel, dPhiR*r, dPhiR*l);
            drawBel(image, bel);        
            printf("done\n");
//            
//            float pointR[3] = {0,0,0};
//            float pointW[3];
//            float npm[3];
//            float newPos[3];
//            if (distF>=0) {
//                
//                robotToSensorPoint(pos, sensorFrontPos, distF, pointW);
//                printf("pointW: [%.2f, %.2f,  %.2f]\n", pointW[0], pointW[1], pointW[2]);
//                plotData(image, pointW);
//                
////                sensorToRobot(distF, sensorFrontPos, pointR);
////                //printf("pointR: [%.2f, %.2f,  %.2f]\n", pointR[0], pointR[1], pointR[2]);
////                robotToWorld(pointR, pos, pointW);
////                //printf("pointW: [%.2f, %.2f,  %.2f]\n", pointW[0], pointW[1], pointW[2]);
////                nearesPointInMap(pointW, npm);
//////                plotLine(image, pointW, npm);
//////                plotData(image, pointW);
//////                robotPosFromWall(distF, npm, pos, newPos);
//////                plotData(image, newPos);
//            }
//            if (distL>=0) {
////                sensorToRobot(distL, sensorLeftPos, pointR);
////                //printf("pointR: [%.2f, %.2f,  %.2f]\n", pointR[0], pointR[1], pointR[2]);
////                robotToWorld(pointR, pos, pointW);
////                //printf("pointW: [%.2f, %.2f,  %.2f]\n", pointW[0], pointW[1], pointW[2]);
////                nearesPointInMap(pointW, npm);
////                plotLine(image, pointW, npm);
////                plotData(image, pointW);
//            }
//            if (distR>=0) {
////                sensorToRobot(distR, sensorRightPos, pointR);
////                //printf("pointR: [%.2f, %.2f,  %.2f]\n", pointR[0], pointR[1], pointR[2]);
////                robotToWorld(pointR, pos, pointW);
////                //printf("pointW: [%.2f, %.2f,  %.2f]\n", pointW[0], pointW[1], pointW[2]);
////                nearesPointInMap(pointW, npm);
////                plotLine(image, pointW, npm);
////                plotData(image, pointW);
//            }
            
            //printf("%.2f | %.2f | %.2f \n", distL, distF, distR);
            
            //Set new target speeds: robot going in a circle:
            simxFloat phiL = 1; //rad/s
            simxFloat phiR = -1; //rad/s            
            setTargetSpeed(clientID, phiL, phiR);

            //Let some time for V-REP do its work:
            extApi_sleepMs(2);
            cv::imshow("Map", image);
            cv::waitKey(25);
        }
        
        //Stop the robot and disconnect from V-Rep;
        //setTargetSpeed(clientID, 0, 0);
        simxStopSimulation(clientID, simx_opmode_oneshot_wait);
        simxFinish(clientID);
        
    } else {
        printf("Nao foi possivel conectar.\n");
        return -2;
    }
    
    return 0;
}


