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

int clientID;
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

#define MAP_SIZE 4.0
#define BEL_NXY 81.0
#define BEL_NTHETA ((360.0/5)+1)
#define XY_STEP (MAP_SIZE/BEL_NXY)
#define TH_STEP (2*M_PI/BEL_NTHETA)

#define IMG_WIDTH 400
#define IMG_BORDER 50

#define SEARCH_SQR 1.0
#define MINPROB 10/(BEL_NXY*BEL_NXY*BEL_NTHETA)

//Robot parameters:
#define r 0.02
#define l 0.1

#define SIGMA_SONAR 0.05
#define SONAR_ANGLE (7.5*2*M_PI/360)
#define SONAR_DELTA (SONAR_ANGLE*2.0/10) 

#define Kr 1
#define Kl 1

float mapLines[20][4];//at most 20 lines
int nLines=0;

//DRAWING MAP
#define DRAW_CLEAR  1
#define DRAW_MAP    1
#define DRAW_BELL   1
#define DRAW_REAL_ROBOT_POS     1
#define DRAW_ESTIM_ROBOT_POS    1
#define DRAW_SONAR_FROM_REAL_POS  1
#define DRAW_SONAR_FROM_ESTIM_POS  1

cv::Mat image(IMG_WIDTH+2*IMG_BORDER, IMG_WIDTH+2*IMG_BORDER, CV_8UC3, cv::Scalar(0,0,0));

cv::Scalar wallColor(255, 0, 0);
cv::Scalar belColor(0,0,255);
cv::Scalar robotColor(0, 255, 0);
cv::Scalar sonarLine(0,255,255);

cv::Scalar lineColor(255,0,255);
cv::Scalar pointColor(255,0,255);

cv::Scalar red(0, 0, 255);
cv::Scalar green(0, 255, 0);
cv::Scalar pink(255, 0, 255);
cv::Scalar yellow(0,255,255);

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

float to180range(float angle) {

   angle = fmod(angle, 2*M_PI);
    if (angle<-M_PI) {
            angle = angle + 2*M_PI;
    }
    else if (angle>M_PI) {
            angle = angle - 2*M_PI;
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

int readOdometers(int clientID, simxFloat &dwL, simxFloat &dwR) {
    static bool first = true;
    //old joint angle position
    static simxFloat lwprev=0; 
    static simxFloat rwprev=0;
    
    //current joint angle position
    simxFloat lwcur=0;
    simxFloat rwcur=0;
    
    if (first) {
        simxInt ret = simxGetJointPosition(clientID, leftMotorHandle, &lwprev, simx_opmode_streaming);
        if (ret>0) return -1;
    
        ret = simxGetJointPosition(clientID, rightMotorHandle, &rwprev, simx_opmode_streaming);
        if (ret>0) return -1;
        
        dwR = dwL = 0;
        first = false;
    }

    simxInt ret = simxGetJointPosition(clientID, leftMotorHandle, &lwcur, simx_opmode_buffer);
    if (ret>0) return -1;
    
    ret = simxGetJointPosition(clientID, rightMotorHandle, &rwcur, simx_opmode_buffer);
    if (ret>0) return -1;
    
    dwL = smallestAngleDiff(lwcur, lwprev);
    dwR = smallestAngleDiff(rwcur, rwprev);
    
//    dwR = rwcur - rwprev;
//    dwL = lwcur - lwprev;
        
    if (fabs(dwR)>M_PI || fabs(dwL)>M_PI) {
        printf("wL: %f - (%f) = %f\n", lwcur, lwprev, dwL);
        printf("wR: %f - (%f) = %f\n", rwcur, rwprev, dwR);
        lwprev = lwcur;
        rwprev = rwcur;
        return -1;
    }
    
    lwprev = lwcur;
    rwprev = rwcur;
    return 0;
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
    
    out[0] = round((vetin[0]+MAP_SIZE/2)*(BEL_NXY-1)/MAP_SIZE);
    out[1] = round((vetin[1]+MAP_SIZE/2)*(BEL_NXY-1)/MAP_SIZE);
    out[2] = round((vetin[2]+M_PI)*(BEL_NTHETA-1)/(2*M_PI));
}

void mapToWorld(float *vetin, float *out) {
    
    out[0] = vetin[0]*MAP_SIZE/(BEL_NXY-1) - MAP_SIZE/2;
    out[1] = vetin[1]*MAP_SIZE/(BEL_NXY-1) - MAP_SIZE/2;
    out[2] = vetin[2]*2*M_PI/(BEL_NTHETA-1) - M_PI;
}

void worldToImage(float *vetin, float *out) {
    
    out[0] = IMG_BORDER+(vetin[0]+(MAP_SIZE/2))*(IMG_WIDTH-1)/MAP_SIZE;
    out[1] = IMG_BORDER+((MAP_SIZE/2) - vetin[1])*(IMG_WIDTH-1)/MAP_SIZE;
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
  
  if (t<0) t = 0;
  if (t>1) t = 1;
  
  npl[0] = pA[0] + AB[0]*t;
  npl[1] = pA[1] + AB[1]*t;
  
//  printf("p: [%.2f, %.2f]\n", p[0], p[1]);
//  printf("pA: [%.2f, %.2f]\n", pA[0], pA[1]);
//  printf("pB: [%.2f, %.2f]\n", pB[0], pB[1]);
//  printf("AP: [%.2f, %.2f]\n", AP[0], AP[1]);
//  printf("AB: [%.2f, %.2f]\n", AB[0], AB[1]);
//  printf("npl: [%.2f, %.2f]\n", npl[0], npl[1]);
//  printf("t: %.2f]\n", t);

  return 1; //point is inside AB
}


//Find nearest point disregarding sonar wave crossing walls
float nearesPointInMap(float *p, float *npl) {
    float pl[2], d, best_d=99999;
    
    for (int i=0; i<nLines; i++) {
        if (nearesPointInLine(p, &mapLines[i][0], &mapLines[i][2], pl))
        {
            d = distance2(pl, p);
            if (d<best_d) {
                best_d = d;
                npl[0] = pl[0];
                npl[1] = pl[1];
            }
//            printf("\np: [%.2f, %.2f] ", p[0], p[1]);
//            printf("pl: [%.2f, %.2f] ", pl[0], pl[1]);
//            printf("d:%.4f bd: %.4f ", d, best_d);
//            printf("npl: [%.2f, %.2f]\n", npl[0], npl[1]);
        }
    }
    
    return best_d;
}

int segmentsIntersect(float *a, float *b, float *c, float *d) {
    
    float den = ((d[1]-c[1])*(b[0]-a[0])-(d[0]-c[0])*(b[1]-a[1]));
    float num1 = ((d[0] - c[0])*(a[1]-c[1]) - (d[1]- c[1])*(a[0]-c[0]));
    float num2 = ((b[0]-a[0])*(a[1]-c[1])-(b[1]-a[1])*(a[0]-c[0]));

    if (den == 0 && num1  == 0 && num2 == 0)
        return -1; // The two lines are coincidents
    if (den == 0)
        return -2; // The two lines are parallel

    float u1 = num1/den;
    float u2 = num2/den;

    if (u1 <0 || u1 > 1 || u2 < 0 || u2 > 1)
        
        return -3; // Lines do not collide 
    
    return 1; // Lines DO collide
}

// Find nearest point considering intersections
// since sonar wave can't go through walls
float nearesPointInMapInter(float *robotPos, float *p, float *npl) {
    float pl[2], d=-1, best_d=99999;
    
    //TODO: check lines de intesect first
    
    for (int i=0; i<nLines; i++) {
        if (nearesPointInLine(p, &mapLines[i][0], &mapLines[i][2], pl))
        {
            d = distance2(pl, p);
            if (d<best_d) {
                best_d = d;
                npl[0] = pl[0];
                npl[1] = pl[1];
            }
//            printf("\np: [%.2f, %.2f] ", p[0], p[1]);
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
    worldToImage(pointW, pointM);
    
    if (pointM[0]>=0 && pointM[0]<image.rows && pointM[1]>=0 && pointM[1]<image.cols) {
        cv::circle(image, cv::Point(pointM[0],pointM[1]), 2, pointColor, 1);
    } else
        printf("Invalid pointM: [%.2f, %.2f]\n", pointM[0], pointM[1]);
}

void plotLine(cv::Mat &image, float* pointWA, float* pointWB) {
    float pointMA[2], pointMB[2];
    worldToImage(pointWA, pointMA);
    worldToImage(pointWB, pointMB);
    
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
        int res;
        while (1) {
            res = fscanf(fp, "%f%f%f%f", &mapLines[nLines][0], &mapLines[nLines][1], &mapLines[nLines][2], &mapLines[nLines][3]);
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

void robotToSensorPoint(float *pRobot, float *pSensor, float dist, float* point){
    float pSR[2];
    sensorToRobot(dist, pSensor, pSR);
    robotToWorld(pSR, pRobot, point);
}

double z(double x) //normal pdf
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

float pGaussian(float dist, float sigma, float step) {
    float x = fabs(dist);
    float halfstep = step/2.0;    
    float zmax = (x+halfstep)/sigma;
    float zmin = (x-halfstep)/sigma;
    return z(zmax) - z(zmin);
}

float pTrvGaussian(cv::Mat &m, cv::Mat &x, cv::Mat E) {
    /*
     * "Approximates" a trivariate gaussian by multiplying tree univariate gaussians with 
     * sigma equal to the main diagonal. It seems to be fast and precise enought for our purposes.
     * 
    */
    return (pGaussian(x.at<float>(0,0)-m.at<float>(0,0), E.at<float>(0,0), XY_STEP)*
            pGaussian(x.at<float>(0,1)-m.at<float>(0,1), E.at<float>(1,1), XY_STEP)*
            pGaussian(smallestAngleDiff(x.at<float>(0,2),m.at<float>(0,2)), E.at<float>(2,2), TH_STEP));///3;
}

void perceptionUpdate(cv::Mat &bel, float distF, float distL, float distR){
    float robotPos[3]; //[x,y,theta]
    int x,y,t;
    float sensorPoint[3], mapPoint[3];
    float p, sum=0;
        
    for(robotPos[0]=-2, x=0; x<BEL_NXY; robotPos[0]+=4.0/(BEL_NXY-1), x++)
    for(robotPos[1]=-2, y=0; y<BEL_NXY; robotPos[1]+=4.0/(BEL_NXY-1), y++)
    for(robotPos[2]=-M_PI, t=0; t<BEL_NTHETA; robotPos[2]+=2*M_PI/(BEL_NTHETA-1), t++) {
        
        float b = bel.at<float>(x,y,t);
        
        if (b>0) {
            float stmin, stmax, d;
            float dF=999, dL=999, dR=999;
            float sensorDir[3];
            
//            robotToSensorPoint(robotPos, sensorFrontPos, distF, sensorPoint);
//            dF = nearesPointInMap(sensorPoint, mapPoint);
            
            sensorDir[0] = sensorFrontPos[0];
            sensorDir[1] = sensorFrontPos[1];
            stmin = sensorFrontPos[2]-SONAR_ANGLE;
            stmax = sensorFrontPos[2]+SONAR_ANGLE;
            for (sensorDir[2] = stmin; sensorDir[2]<=stmax; sensorDir[2]+=SONAR_DELTA) {
                robotToSensorPoint(robotPos, sensorDir, distF, sensorPoint);
                d = nearesPointInMap(sensorPoint, mapPoint);
                if (d<dF) dF=d;
            }
            
//            robotToSensorPoint(robotPos, sensorLeftPos, distL, sensorPoint);
//            dL = nearesPointInMap(sensorPoint, mapPoint);
            sensorDir[0] = sensorLeftPos[0];
            sensorDir[1] = sensorLeftPos[1];
            stmin = sensorLeftPos[2]-SONAR_ANGLE;
            stmax = sensorLeftPos[2]+SONAR_ANGLE;
            for (sensorDir[2] = stmin; sensorDir[2]<=stmax; sensorDir[2]+=SONAR_DELTA) {
                robotToSensorPoint(robotPos, sensorDir, distL, sensorPoint);
                d = nearesPointInMap(sensorPoint, mapPoint);
                if (d<dL) dL=d;
            }

//            robotToSensorPoint(robotPos, sensorRightPos, distR, sensorPoint);
//            dR = nearesPointInMap(sensorPoint, mapPoint);
            sensorDir[0] = sensorRightPos[0];
            sensorDir[1] = sensorRightPos[1];
            stmin = sensorRightPos[2]-SONAR_ANGLE;
            stmax = sensorRightPos[2]+SONAR_ANGLE;
            for (sensorDir[2] = stmin; sensorDir[2]<=stmax; sensorDir[2]+=SONAR_DELTA) {
                robotToSensorPoint(robotPos, sensorDir, distR, sensorPoint);
                d = nearesPointInMap(sensorPoint, mapPoint);
                if (d<dR) dR=d;
            }

            p = pGaussian(dF, SIGMA_SONAR, XY_STEP)*pGaussian(dL, SIGMA_SONAR, XY_STEP)*pGaussian(dR, SIGMA_SONAR, XY_STEP)*b;
            bel.at<float>(x,y,t) = p;
            sum+=p;
        }
    }
    
    if (sum>0)
        bel = bel/sum;
//    printf("robot bel(%d,%d,%d): %f\n", xr,yr,tr, bel.at<float>(xr,yr,tr));
//    printf("max p bel(%d,%d,%d): %f\n", xp,yp,tp, bel.at<float>(xp,yp,tp));
}


void fdeltaRL(float theta, float ds, float dtheta, cv::Mat &FDrl) {
    //ds = fabs(ds);
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
    Ep = FDrl*Ed*FDrl.t();

//    cout << "Ed:\n" << Ed.at<float>(0,0) << " | "  << Ed.at<float>(0,1) << endl;
//    cout << Ed.at<float>(1,0) << " | "  << Ed.at<float>(1,1) << endl;
//    
//    cout << "FDrl:\n" << FDrl.at<float>(0,0) << " | "  << FDrl.at<float>(0,1) << endl;
//    cout << FDrl.at<float>(1,0) << " | "  << FDrl.at<float>(1,1) << endl;
//    cout << FDrl.at<float>(2,0) << " | "  << FDrl.at<float>(2,1) << endl;
//    
//    cout << "Ep:\n" << Ep.at<float>(0,0) << " | "  << Ep.at<float>(0,1) << " | "  << Ep.at<float>(0,2) << endl;
//    cout << Ep.at<float>(1,0) << " | "  << Ep.at<float>(1,1) << " | "  << Ep.at<float>(1,2) << endl;
//    cout << Ep.at<float>(2,0) << " | "  << Ep.at<float>(2,1) << " | "  << Ep.at<float>(2,2) << endl; 
}

int actionUpdate(cv::Mat &bel, float dsl, float dsr) {
    static float old_dsl = 0;
    static float old_dsr = 0;

    //printf("\ndsr: %f, dsl: %f\n", dsr, dsl);
    dsr = old_dsr = (old_dsr+dsr);
    dsl = old_dsl = (old_dsl+dsl);
    float ds = (dsr + dsl)/2;
    float dtheta = (dsr - dsl)/(2*l);
    //printf("ndsr: %f, ndsl: %f, abs ds:%f<%f, dt: %f<%f\n", dsr, dsl, fabs(ds), XY_STEP, fabs(dtheta), TH_STEP);
    
    if (fabs(ds) < XY_STEP && fabs(dtheta) < TH_STEP) {
        return -1;
    } else {
        old_dsl = 0;
        old_dsr = 0;        
    }
    
    int x1,y1,t1; //x,y,theta indexes on bel map
    int x0,y0,t0; //x,y,theta indexes on tempbel map
    
    cv::Mat X1(1,3, CV_32FC1), //current position [x,y,theta]
            X0(1,3, CV_32FC1), //previous position [x,y,theta]
            M(1,3, CV_32FC1),  //new expected position [x,y,theta]
            Ed(2,2, CV_32FC1),  //covar(dsr, dsl)
            FDrl(3,2, CV_32FC1),//Jacobian
            Ep(3,3, CV_32FC1);  //motion error Sigmap
    
    const int belSizes[3]={BEL_NXY, BEL_NXY, BEL_NTHETA};
    cv::Mat sumbel(3, belSizes, CV_32FC1, 0.0); //somatory belief mat init with zeroes
        
    Ed.at<float>(0,0) = Kr*fabs(dsr);
    Ed.at<float>(1,1) = Kl*fabs(dsl);
    Ed.at<float>(0,1) = 0;
    Ed.at<float>(1,0) = 0;
    
    int cells=0;
    float sum=0;
    float costdt2=0;
    float sintdt2=0;
    
    for(X0.at<float>(0,2)=-M_PI, t0=0; t0<BEL_NTHETA; X0.at<float>(0,2)+=2*M_PI/(BEL_NTHETA-1), t0++) {//for each old theta
        bool newTheta = true;
        
        for(X0.at<float>(0,0)=-2, x0=0; x0<BEL_NXY; X0.at<float>(0,0)+=4.0/(BEL_NXY-1), x0++) { //for each old x
            for(X0.at<float>(0,1)=-2, y0=0; y0<BEL_NXY; X0.at<float>(0,1)+=4.0/(BEL_NXY-1), y0++) {  //for each old y

                float b = bel.at<float>(x0,y0,t0);
                if (b>MINPROB) {
                    
                    if (newTheta) {
                        costdt2 = cos(X0.at<float>(0,2)+dtheta/2);
                        sintdt2 = sin(X0.at<float>(0,2)+dtheta/2);
                        fdeltaRL(X0.at<float>(0,2), ds, dtheta, FDrl);
                        odomError(FDrl, Ed, Ep);
                        
                        M.at<float>(0,2) = to180range(X0.at<float>(0,2)+dtheta);
                        newTheta = false;
                    }
                    
                    M.at<float>(0,1) = X0.at<float>(0,1)+ds*sintdt2;
                    M.at<float>(0,0) = X0.at<float>(0,0)+ds*costdt2;
                    
                    //Update only a rectangle around X0
                    float dsrect = SEARCH_SQR*(fabs(ds)+XY_STEP);
                    float dtrect = SEARCH_SQR*(fabs(dtheta)+TH_STEP);//2*(dtheta + 1 step) since dtheta can be 0
                    float x1min[3], x1max[3];
                    //x1min
                    x1min[0] = M.at<float>(0,0)-dsrect;
                    if (x1min[0]<-2) x1min[0] = -2;                    
                    x1min[1] = M.at<float>(0,1)-dsrect;
                    if (x1min[1]<-2) x1min[1] = -2;
                    x1min[2] = M.at<float>(0,2)-dtrect;
                    if (x1min[2]<-M_PI) x1min[2] = -M_PI;
                    
                    //x1max
                    x1max[0] = M.at<float>(0,0)+dsrect;
                    if (x1max[0]>2) x1max[0] = 2;                    
                    x1max[1] = M.at<float>(0,1)+dsrect;
                    if (x1max[1]>2) x1max[1] = 2;
                    x1max[2] = M.at<float>(0,2)+dtrect;
                    if (x1max[2]>M_PI) x1max[2] = M_PI;
                    
                    //get map coordinates
                    float mapmin[3], mapmax[3];
                    worldToMap(x1min, mapmin);
                    worldToMap(x1max, mapmax);     
                    
                    mapToWorld(mapmin, x1min);
                    mapToWorld(mapmax, x1max);
                    
                    
//                    printf("X0 = (%.2f,%.2f,%.2f) M = (%.2f,%.2f,%.2f)\n", X0.at<float>(0,0), X0.at<float>(0,1), X0.at<float>(0,2), M.at<float>(0,0), M.at<float>(0,1), M.at<float>(0,2));
                    
                    for(X1.at<float>(0,2) = x1min[2], t1 =  mapmin[2]; t1<=mapmax[2]; X1.at<float>(0,2)+=2*M_PI/(BEL_NTHETA-1), t1++) {//for each new theta
                        for(X1.at<float>(0,0) = x1min[0], x1 =  mapmin[0]; x1<=mapmax[0]; X1.at<float>(0,0)+=4.0/(BEL_NXY-1), x1++) { //for each new x        
                            for(X1.at<float>(0,1) = x1min[1], y1 =  mapmin[1]; y1<=mapmax[1]; X1.at<float>(0,1)+=4.0/(BEL_NXY-1), y1++) { //for each new y           

                                float px1_u1x0 = pTrvGaussian(M, X1, Ep);
                                sumbel.at<float>(x1,y1,t1) += b*px1_u1x0;
                                sum += b*px1_u1x0;
                                cells++;
//                                if (px1_u1x0>0) {
//                                    printf("px1_u1x0(%.2f,%.2f,%.2f | %.2f,%.2f,%.2f): p:%f  b:%f pb:%f\n", X1.at<float>(0,0), X1.at<float>(0,1), X1.at<float>(0,2), M.at<float>(0,0), M.at<float>(0,1), M.at<float>(0,2), px1_u1x0, b, px1_u1x0*b);
//                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    if (sum>0) {
        bel = sumbel/sum;   
    }
    else { 
        printf("I'm lost!!!\n");
        bel = 1.0; //start over with perception data only
//        sum = cv::sum(bel)[0];
//        bel = bel/sum;
        return BEL_NXY*BEL_NXY*BEL_NTHETA;
    }
  
//    double testMaxval;
//    int maxIdx[3];
//    cv::minMaxIdx(sumbel, 0, &testMaxval, 0, maxIdx);
//    printf("Max bel: %f (%d,%d,%d) Sum: %f\n", testMaxval, maxIdx[0], maxIdx[1], maxIdx[2], sum); 
//    
    return cells;
}

int drawMap(cv::Mat &image) {

    lineColor = wallColor;
    for (int i=0; i<nLines; i++) {
        plotLine(image, &mapLines[i][0], &mapLines[i][2]);
    }
}

int drawSonar(cv::Mat &image, float robotPos[3], float distF, float distL, float distR) {
    float stmin, stmax, d;
    float dF=999, dL=999, dR=999;
    float sensorDir[3];
    float sensorPoint[3], mapPoint[3], bestPoint[3];
    lineColor = sonarLine;
    
    sensorDir[0] = sensorFrontPos[0];
    sensorDir[1] = sensorFrontPos[1];
    stmin = sensorFrontPos[2]-SONAR_ANGLE;
    stmax = sensorFrontPos[2]+SONAR_ANGLE;
    for (sensorDir[2] = stmin; sensorDir[2]<=stmax; sensorDir[2]+=SONAR_DELTA) {
        robotToSensorPoint(robotPos, sensorDir, distF, sensorPoint);
        d = nearesPointInMap(sensorPoint, mapPoint);
        if (d<dF) {
            dF=d;
            bestPoint[0] = mapPoint[0];
            bestPoint[1] = mapPoint[1];
            bestPoint[2] = mapPoint[2];
        }
    }
    
    plotLine(image, robotPos, mapPoint);

    sensorDir[0] = sensorLeftPos[0];
    sensorDir[1] = sensorLeftPos[1];
    stmin = sensorLeftPos[2]-SONAR_ANGLE;
    stmax = sensorLeftPos[2]+SONAR_ANGLE;
    for (sensorDir[2] = stmin; sensorDir[2]<=stmax; sensorDir[2]+=SONAR_DELTA) {
        robotToSensorPoint(robotPos, sensorDir, distL, sensorPoint);
        d = nearesPointInMap(sensorPoint, mapPoint);
        if (d<dL) {
            dL=d;
            bestPoint[0] = mapPoint[0];
            bestPoint[1] = mapPoint[1];
            bestPoint[2] = mapPoint[2];
        }
    }
    
    plotLine(image, robotPos, mapPoint);

    sensorDir[0] = sensorRightPos[0];
    sensorDir[1] = sensorRightPos[1];
    stmin = sensorRightPos[2]-SONAR_ANGLE;
    stmax = sensorRightPos[2]+SONAR_ANGLE;
    for (sensorDir[2] = stmin; sensorDir[2]<=stmax; sensorDir[2]+=SONAR_DELTA) {
        robotToSensorPoint(robotPos, sensorDir, distR, sensorPoint);
        d = nearesPointInMap(sensorPoint, mapPoint);
        if (d<dR)  {
            dR=d;
            bestPoint[0] = mapPoint[0];
            bestPoint[1] = mapPoint[1];
            bestPoint[2] = mapPoint[2];
        }
    }
    
    plotLine(image, robotPos, mapPoint);
}

void drawBel(cv::Mat &image, cv::Mat bel) {
    float robotPos[3]; //[x,y,theta]
    int x,y,t;
    float p, min=9999, max=0;
    
    for(x=0; x<BEL_NXY; x++)
    for(y=0; y<BEL_NXY; y++)
    for(t=0; t<BEL_NTHETA;  t++) {
        p=bel.at<float>(x, y, t);
        if (p>max) max = p;
        if (p<min) min = p;
    }
            
    for(robotPos[0]=-2, x=0; x<BEL_NXY; robotPos[0]+=4.0/(BEL_NXY-1), x++)
    for(robotPos[1]=-2, y=0; y<BEL_NXY; robotPos[1]+=4.0/(BEL_NXY-1), y++) {
        float p=0;
        for(t=0; t<BEL_NTHETA;  t++) {
            if (bel.at<float>(x, y, t)>p)
                p = bel.at<float>(x, y, t);
        }

        int color = (p-min)/(max-min)*255;      
        
        pointColor = cv::Scalar(color, color, color);
        plotData(image, robotPos);
    }
}

void drawRobot(cv::Mat &image, float robotPos[3]) {
    float pointM[2];
    worldToImage(robotPos, pointM);
    
    cv::circle(image, cv::Point(pointM[0],pointM[1]), 0.7*BEL_NXY/4, robotColor, 2);
}

void draw(cv::Mat &bel, float distF, float distL, float distR) {
    
    if (DRAW_CLEAR)
        image = 0.0; //clear image
    
    if (DRAW_MAP)
        drawMap(image);
    
    if (DRAW_BELL)
        drawBel(image, bel);
    
    if (DRAW_REAL_ROBOT_POS) {
        float pos[3];
        getPosition(clientID, pos);
        
        robotColor = green;
        drawRobot(image, pos);
        
        sonarLine = pink;        
        if (DRAW_SONAR_FROM_REAL_POS && distF>=0 && distL>=0 && distR>=0) {
            drawSonar(image, pos, distF, distL, distR);
        }
    }
    
    if (DRAW_ESTIM_ROBOT_POS) {
        int maxIdx[3];
        double testMaxval;
        float pos[3];
        cv::minMaxIdx(bel, 0, &testMaxval, 0, maxIdx);
        float idxfloat[3] = {maxIdx[0], maxIdx[1], maxIdx[2]};
        mapToWorld(idxfloat, pos);
        
        robotColor = red;
        drawRobot(image, pos);
        
        sonarLine = yellow;
        if (DRAW_SONAR_FROM_ESTIM_POS && distF>=0 && distL>=0 && distR>=0) {
            drawSonar(image, pos, distF, distL, distR);
        }
    }
    
    cv::imshow("Map", image);
    cv::waitKey(25);
}

int main(int argc, char* argv[]) {

    std::string ipAddr = V_REP_IP_ADDRESS;
    int portNb = V_REP_PORT;

    if (argc > 1) {
        ipAddr = argv[1];
    }

    printf("Iniciando conexao com: %s...\n", ipAddr.c_str());

    clientID = simxStart((simxChar*) (simxChar*) ipAddr.c_str(), portNb, true, true, 2000, 5);
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

        if (!loadMap("map.txt")) {
            printf("Não foi possível carregar o mapa.\n");
            return -1;
        }
        
        printf("Simulação iniciada.\n");

        simxFloat detectedPoint[3];
        simxUChar detectionState=1;
        simxReadProximitySensor(clientID, sensorFrontHandle, &detectionState, detectedPoint, NULL, NULL, simx_opmode_streaming);
        simxReadProximitySensor(clientID, sensorLeftHandle, &detectionState, detectedPoint, NULL, NULL, simx_opmode_streaming);
        simxReadProximitySensor(clientID, sensorRightHandle, &detectionState, detectedPoint, NULL, NULL, simx_opmode_streaming);
        
        const int belSizes[3]={BEL_NXY, BEL_NXY, BEL_NTHETA};
        cv::Mat bel(3, belSizes, CV_32FC1);

       /*
        * Initialize bel with 100% certainty at real robot position so
        * that the tracking starts correctly.
        */ 
//        Get real robot position from v-rep API:
        float pos[3] = {0,0,0};
        getPosition(clientID, pos); //comment out to simulate robot in position {0,0,0} for debug        
        float firstMapPos[3];
        worldToMap(pos, firstMapPos);
        bel = 0.0;
        bel.at<float>(firstMapPos[0],firstMapPos[1], firstMapPos[2]) = 1;

        //start simulation
        int ret = simxStartSimulation(clientID, simx_opmode_oneshot_wait);
        
        if (ret==-1) {
            printf("Não foi possível iniciar a simulação.\n");
            return -1;
        }
        
        int cells=1;
        //While is connected:
        while (simxGetConnectionId(clientID) != -1)
        {  
             // Perception update //
            simxFloat distF=-1, distR=-1, distL=-1;
            distF = readSonar(clientID, sensorFrontHandle);
            distL = readSonar(clientID, sensorLeftHandle);
            distR = readSonar(clientID, sensorRightHandle);

            if (cells>0 && distF>=0 && distL>=0 && distR>=0) {
                //Do only if action step was successfull and sonar data is good
                //printf("start perception update...");fflush(stdout);
                perceptionUpdate(bel, distF, distL, distR);
                printf("perception update done\n");
            }

            // Action update //
            //Read current wheels angle variation:
            simxFloat dwL, dwR; //rad
            if (readOdometers(clientID, dwL, dwR)==0) {
                //printf("start action update width: %.4f %.4f...", dwL*r, dwR*r);fflush(stdout);
                cells = actionUpdate(bel, dwL*r, dwR*r);
                if (cells>0)
                    printf("action update %d cells\n", cells);
            } else
                printf("Error reading odometers\n");

            //Set robot target speed if you wish:
//            simxFloat phiL = 0;//1; //rad/s
//            simxFloat phiR = 0;//-1;//rad/s            
//            setTargetSpeed(clientID, phiL, phiR);

            //draw our debug diagram:
            draw(bel, distF, distL, distR);

            //Read simulation time of the last command:
            simxInt time = getSimTimeMs(clientID); //Simulation time in ms or 0 if sim is not running
            //stop the loop if simulation is has been stopped:
            if (time == 0) break;                
            
            //Let some time for V-REP do its work:
            extApi_sleepMs(2);            
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