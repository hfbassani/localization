#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ffns.h"

double d_nint(doublereal *x)
{
return( (*x)>=0 ?
	floor(*x + .5) : -floor(.5 - *x) );
}

double d_int(doublereal *x)
{
return( (*x>0) ? floor(*x) : -floor(- *x) );
}

double pow_di(doublereal *ap, integer *bp)
{
double pow, x;
integer n;

pow = 1;
x = *ap;
n = *bp;

if(n != 0)
	{
	if(n < 0)
		{
		n = -n;
		x = 1/x;
		}
	for( ; ; )
		{
		if(n & 01)
			pow *= x;
		if(n >>= 1)
			x *= x;
		else
			break;
		}
	}
return(pow);
}

double d_sign(doublereal *a, doublereal *b)
{
double x;
x = (*a >= 0 ? *a : - *a);
return( *b >= 0 ? x : -x);
}


double pow_dd(doublereal *ap, doublereal *bp)
{
return(pow(*ap, *bp) );
}

double d_lg10(doublereal *x)
{
return( log10e * log(*x) );
}

/*************** Common Block Declarations */
static union {
    struct {
	doublereal b2p, b3p, rho21, rho31, rho;
    } _1;
    struct {
	doublereal b2, b3, rho21, rho31, rho;
    } _2;
} trvbkd_;

#define trvbkd_1 (trvbkd_._1)
#define trvbkd_2 (trvbkd_._2)

/* Table of constant values */
static doublereal c_b5 = -8.5;
static doublereal c_b6 = 5e-16;

doublereal bvn_(real *sh, real *sk, real *r__)
{
    /* Initialized data */
    static real x[5] = { (float)-.9061798,(float)-.5384693,(float)0. };
    static real w[5] = { (float).2369269,(float).4786287,(float).5688889 };

    /* System generated locals */
    real ret_val, r__1, r__2, r__3, r__4;

    /* Local variables */
    static real a, b, c__, h__;
    static integer i__;
    static real k, as, bs, hk, hs, sn, rs, xs;
    static real asr;


/*     A function for computing bivariate normal probabilities. */
/*     This function uses an algorithm given in the paper */
/*        "Numerical Computation of Bivariate and */
/*         Trivariate Normal Probabilities", by: */
/*       Alan Genz */
/*       Department of Mathematics */
/*       Washington State University */
/*       Pullman, WA 99164-3113 */
/*       Email : alangenz@wsu.edu */

/* BVN - calculate the probability that X is larger than SH and Y is */
/*       larger than SK. */

/* Parameters */
/*   SH  REAL, integration limit */
/*   SK  REAL, integration limit */
/*   R   REAL, correlation coefficient */
/*   LG  INTEGER, number of Gauss Rule Points and Weights */

/*   X(I) = +-SQRT((5+-SQRT(10/7))/9), 0 */
    x[3] = -x[1];
    x[4] = -x[0];
    w[3] = w[1];
    w[4] = w[0];
    h__ = *sh;
    k = *sk;
    hk = h__ * k;
    ret_val = (float)0.;

/*     Absolute value of the correlation coefficient is less than 0.8 */
    if (dabs(*r__) < (float).8) {
	hs = (h__ * h__ + k * k) / 2;
	asr = (float)asin(*r__) / 2;
	for (i__ = 1; i__ <= 5; ++i__) {
	    sn = (float)sin(asr * (x[i__ - 1] + 1));
	    ret_val += (float)(w[i__ - 1] * exp((sn * hk - hs) / (1 - sn * sn)));
	}
	r__1 = -h__;
	r__2 = -k;
	ret_val = (float)(asr * ret_val / (float)6.283185 + phi_(&r__1) * phi_(&r__2));
    } else {

/*     For larger correlation coefficient */
	if (*r__ < (float)0.) {
	    k = -k;
	    hk = -hk;
	}
	if (dabs(*r__) < (float)1.) {
	    as = (1 - *r__) * (*r__ + 1);
	    a = (float)sqrt(as);
/* Computing 2nd power */
	    r__1 = h__ - k;
	    bs = r__1 * r__1;
	    c__ = (4 - hk) / 8;
	    ret_val = (float)(a * exp(-(bs / as + hk) / 2) * (1 - c__ * (bs - as) / 3));
	    if (-hk < (float)100.) {
		   b = (float)sqrt(bs);
		   r__1 = -b / a;
		   ret_val -= (float)(exp(-hk / 2) * sqrt((float)6.283185) * phi_(&r__1) 
			* b * (1 - c__ * bs / 3));}
	    a /= 2;
	    for (i__ = 1; i__ <= 5; ++i__) {
/* Computing 2nd power */
		r__1 = a * (x[i__ - 1] + 1);
		xs = r__1 * r__1;
		rs = (float)sqrt(1 - xs);
		asr = -(bs / xs + hk) / 2;
		if (asr > (float)-100.) {
		    ret_val += (float)(a * w[i__ - 1] * (exp(-bs / (xs * 2) - hk / (
			    rs + 1)) / rs - exp(asr) * (c__ * xs + 1)));}
	    }
	    ret_val = -ret_val / (float)6.283185;
	}
	if (*r__ > (float)0.) {
	    r__1 = (float)(-dmax(h__,k));
	    ret_val += (float)phi_(&r__1);
	}
	if (*r__ < (float)0.) {
/* Computing MAX */
	    r__3 = -h__;
	    r__4 = -k;
	    r__1 = (float)0., r__2 = (float)(phi_(&r__3) - phi_(&r__4));
	    ret_val = (float)(-ret_val + dmax(r__1,r__2));}
    }
    return ret_val;
} /* bvn_ */


doublereal phi_(real *x)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1;

    /* Local variables */
    d__1 = (doublereal) (*x);
    ret_val = (float)phid_(&d__1);
    return ret_val;
} /* phi_ */


doublereal bvnd_(doublereal *dh, doublereal *dk, doublereal *r__)
{
    /* Initialized data */
    static struct {
	doublereal e_1[3];
	doublereal fill_2[7];
	doublereal e_3[6];
	doublereal fill_4[4];
	doublereal e_5[10];
	} equiv_38 = { .1713244923791705, .3607615730481384, 
		.4679139345726904, {0}, .04717533638651177, .1069393259953183,
		 .1600783285433464, .2031674267230659, .2334925365383547, 
		.2491470458134029, {0}, .01761400713915212, 
		.04060142980038694, .06267204833410906, .08327674157670475, 
		.1019301198172404, .1181945319615184, .1316886384491766, 
		.1420961093183821, .1491729864726037, .1527533871307259 };

#define w ((doublereal *)&equiv_38)

    static struct {
	doublereal e_1[3];
	doublereal fill_2[7];
	doublereal e_3[6];
	doublereal fill_4[4];
	doublereal e_5[10];
	} equiv_39 = { -.9324695142031522, -.6612093864662647, 
		-.238619186083197, {0}, -.9815606342467191, -.904117256370475,
		 -.769902674194305, -.5873179542866171, -.3678314989981802, 
		-.1252334085114692, {0}, -.9931285991850949, 
		-.9639719272779138, -.9122344282513259, -.8391169718222188, 
		-.7463319064601508, -.636053680726515, -.5108670019508271, 
		-.3737060887154196, -.2277858511416451, -.07652652113349733 };

#define x ((doublereal *)&equiv_39)


    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    static doublereal a, b, c__, d__, h__;
    static integer i__;
    static doublereal k;
    static integer lg;
    static doublereal as;
    static integer ng;
    static doublereal bs, hk, hs;
    static integer is;
    static doublereal sn, rs, xs, bvn, asr;

/*     A function for computing bivariate normal probabilities. */
/*       Alan Genz */
/*       Department of Mathematics */
/*       Washington State University */
/*       Pullman, WA 99164-3113 */
/*       Email : alangenz@wsu.edu */

/*    This function is based on the method described by */
/*        Drezner, Z and G.O. Wesolowsky, (1989), */
/*        On the computation of the bivariate normal inegral, */
/*        Journal of Statist. Comput. Simul. 35, pp. 101-107, */
/*   with major modifications for double precision, and for |R| close to 1 */
  
/* BVND - calculate the probability that X is larger than DH and Y is */
/*       larger than DK. */

/* Parameters */
/*   DH  DOUBLE PRECISION, integration limit */
/*   DK  DOUBLE PRECISION, integration limit */
/*   R   DOUBLE PRECISION, correlation coefficient */

/*     Gauss Legendre Points and Weights, N =  6 */
/*     Gauss Legendre Points and Weights, N = 12 */
/*     Gauss Legendre Points and Weights, N = 20 */
    if (abs(*r__) < .3) {
		ng = 1;
		lg = 3;}
		else if (abs(*r__) < .75) {
			ng = 2;
			lg = 6;}
				else {
				ng = 3;
				lg = 10;}
    h__ = *dh;
    k = *dk;
    hk = h__ * k;
    bvn = 0.;
    if (abs(*r__) < .925) {
	hs = (h__ * h__ + k * k) / 2;
	asr = asin(*r__);
	i__1 = lg;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    for (is = -1; is <= 1; is += 2) {
		sn = sin(asr * (is * x[i__ + ng * 10 - 11] + 1) / 2);
		bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * 
			sn));}
	}
	d__1 = -h__;
	d__2 = -k;
	bvn = bvn * asr / 12.566370614359172 + phid_(&d__1) * phid_(&d__2);
    } else {
	if (*r__ < 0.) {
	    k = -k;
	    hk = -hk;
	}
	if (abs(*r__) < 1.) {
	    as = (1 - *r__) * (*r__ + 1);
	    a = sqrt(as);
/* Computing 2nd power */
	    d__1 = h__ - k;
	    bs = d__1 * d__1;
	    c__ = (4 - hk) / 8;
	    d__ = (12 - hk) / 16;
	    asr = -(bs / as + hk) / 2;
	    if (asr > -100.) {
		bvn = a * exp(asr) * (1 - c__ * (bs - as) * (1 - d__ * bs / 5)
			 / 3 + c__ * d__ * as * as / 5);}
	    if (-hk < 100.) {
		b = sqrt(bs);
		d__1 = -b / a;
		bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * phid_(&d__1) *
			 b * (1 - c__ * bs * (1 - d__ * bs / 5) / 3);}
	    a /= 2;
	    i__1 = lg;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		for (is = -1; is <= 1; is += 2) {
/* Computing 2nd power */
		    d__1 = a * (is * x[i__ + ng * 10 - 11] + 1);
		    xs = d__1 * d__1;
		    rs = sqrt(1 - xs);
		    asr = -(bs / xs + hk) / 2;
		    if (asr > -100.) {
			bvn += a * w[i__ + ng * 10 - 11] * exp(asr) * (exp(
				-hk * (1 - rs) / ((rs + 1) * 2)) / rs - (c__ *
				 xs * (d__ * xs + 1) + 1));}
		}
	    }
	    bvn = -bvn / 6.283185307179586;
	}
	if (*r__ > 0.) {
	    d__1 = -max(h__,k);
	    bvn += phid_(&d__1);
	}
	if (*r__ < 0.) {
/* Computing MAX */
	    d__3 = -h__;
	    d__4 = -k;
	    d__1 = 0., d__2 = phid_(&d__3) - phid_(&d__4);
	    bvn = -bvn + max(d__1,d__2);
	}
    }
    ret_val = bvn;
    return ret_val;
} /* bvnd_ */

#undef x
#undef w



doublereal tvnd_(doublereal *limit, doublereal *sigma)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static logical tail;
    static doublereal rho32, b1, b2, b3;
    static doublereal sq21, sq31;

/*     A function for computing trivariate normal probabilities. */
/*     This function uses an algorithm given in the paper */
/*        "Numerical Computation of Bivariate and */
/*             Trivariate Normal Probabilities", by : */
/*       Alan Genz */
/*       Department of Mathematics */
/*       Washington State University */
/*       Pullman, WA 99164-3113 */
/*       Email : alangenz@wsu.edu */

/* TVND calculates the probability that X(I) < LIMIT(I), I = 1, 2, 3. */

/* Parameters */
/*   LIMIT  DOUBLE PRECISION array of three upper integration limits. */
/*   SIGMA  DOUBLE PRECISION array of three correlation coefficients, */
/*          SIGMA should contain the lower left portion of the */
/*          correlation matrix R. */
/*          SIGMA(1) = R(2,1), SIGMA(2) = R(3,1), SIGMA(3) = R(3,2). */

/*    TVND cuts the outer integral over -infinity to B1 to */
/*      an integral from -8.5 to B1 and then uses an adaptive */
/*      integration method to compute the integral of a bivariate */
/*      normal distribution function. */

/*     Bivariate normal distribution function BVND is required. */

    /* Parameter adjustments */
    --sigma;
    --limit;

    /* Function Body */
    b1 = limit[1];
    b2 = limit[2];
    b3 = limit[3];
    trvbkd_1.rho21 = sigma[1];
    trvbkd_1.rho31 = sigma[2];
    rho32 = sigma[3];
/* Computing MAX */
    d__1 = abs(b1), d__2 = abs(b3);
    if (abs(b2) >= max(d__1,d__2)) {
	b1 = b2;
	b2 = limit[1];
	trvbkd_1.rho31 = rho32;
	rho32 = sigma[2];
    } else /* if(complicated condition) */ {
/* Computing MAX */
	d__1 = abs(b1), d__2 = abs(b2);
	if (abs(b3) >= max(d__1,d__2)) {
	    b1 = b3;
	    b3 = limit[1];
	    trvbkd_1.rho21 = rho32;
	    rho32 = sigma[1];
	}
    }
    tail = FALSE_;
    if (b1 > 0.) {
	tail = TRUE_;
	b1 = -b1;
	trvbkd_1.rho21 = -trvbkd_1.rho21;
	trvbkd_1.rho31 = -trvbkd_1.rho31;
    }
    if (b1 > -8.5) {
	if (abs(trvbkd_1.rho21) * 2 < 1.) {
/* Computing 2nd power */
	    d__1 = trvbkd_1.rho21;
	    sq21 = sqrt(1 - d__1 * d__1);
	} else {
	    sq21 = sqrt((1 - trvbkd_1.rho21) * (trvbkd_1.rho21 + 1));
	}
	if (abs(trvbkd_1.rho31) * 2 < 1.) {
/* Computing 2nd power */
	    d__1 = trvbkd_1.rho31;
	    sq31 = sqrt(1 - d__1 * d__1);
	} else {
	    sq31 = sqrt((1 - trvbkd_1.rho31) * (trvbkd_1.rho31 + 1));
	}
	trvbkd_1.rho = (rho32 - trvbkd_1.rho21 * trvbkd_1.rho31) / (sq21 * 
		sq31);
	trvbkd_1.b2p = b2 / sq21;
	trvbkd_1.rho21 /= sq21;
	trvbkd_1.b3p = b3 / sq31;
	trvbkd_1.rho31 /= sq31;
	ret_val = adoned_(&c_b5, &b1, &c_b6) / 2.506628274631;
    } else {
	ret_val = 0.;
    }
    if (tail) {
	d__1 = -b2;
	d__2 = -b3;
	ret_val = bvnd_(&d__1, &d__2, &rho32) - ret_val;
    }
    return ret_val;
} /* tvnd_ */


doublereal trvfnd_(doublereal *t)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    d__1 = *t * trvbkd_2.rho21 - trvbkd_2.b2;
    d__2 = *t * trvbkd_2.rho31 - trvbkd_2.b3;
    ret_val = exp(-(*t) * *t / 2) * bvnd_(&d__1, &d__2, &trvbkd_2.rho);
    return ret_val;
} /* trvfnd_ */


doublereal adoned_(doublereal *a, doublereal *b, doublereal *tol)
//trvfnd_
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal ai[100], bi[100], ei[100], fi[100];
    static integer im, ip;
    static doublereal fin, err;


/*     One Dimensional Globally Adaptive Integration Function */
    ip = 1;
    ai[ip - 1] = *a;
    bi[ip - 1] = *b;
    fi[ip - 1] = krnrdd_(&ai[ip - 1], &bi[ip - 1], &ei[ip - 1]);
    im = 1;
L10:
    ++im;
    bi[im - 1] = bi[ip - 1];
    ai[im - 1] = (ai[ip - 1] + bi[ip - 1]) / 2;
    bi[ip - 1] = ai[im - 1];
    fin = fi[ip - 1];
    fi[ip - 1] = krnrdd_(&ai[ip - 1], &bi[ip - 1], &ei[ip - 1]);
    fi[im - 1] = krnrdd_(&ai[im - 1], &bi[im - 1], &ei[im - 1]);
    err = (d__1 = fin - fi[ip - 1] - fi[im - 1], abs(d__1)) / 2;
    ei[ip - 1] += err;
    ei[im - 1] += err;
    ip = 1;
    err = 0.;
    fin = 0.;
    i__1 = im;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ei[i__ - 1] > ei[ip - 1]) {
	    ip = i__;
	}
	fin += fi[i__ - 1];
	err += ei[i__ - 1];
    }
    if (err > *tol && im < 100) {
	goto L10;
    }
    ret_val = fin;
    return ret_val;
} /* adoned_ */



doublereal krnrdd_(doublereal *a, doublereal *b, doublereal *abserr)
{
    /* Initialized data */
    static doublereal wg[7] = { .2729250867779007,.05566856711617449,
	    .1255803694649048,.1862902109277352,.2331937645919914,
	    .2628045445102478 };
    static doublereal xgk[12] = { 0.,.9963696138895427,.978228658146057,
	    .9416771085780681,.8870625997680953,.8160574566562211,
	    .7301520055740492,.6305995201619651,.5190961292068118,
	    .3979441409523776,.269543155952345,.1361130007993617 };
    static doublereal wgk[12] = { .1365777947111183,.00976544104596129,
	    .02715655468210443,.04582937856442671,.06309742475037484,
	    .07866457193222764,.09295309859690074,.1058720744813894,
	    .1167395024610472,.1251587991003195,.1312806842298057,
	    .1351935727998845 };

    /* System generated locals */
    doublereal ret_val, d__1, d__2;

    /* Local variables */
    static integer j;
    static doublereal fc, abscis, hflgth, center, resltg, resltk, funsum;


/*     Kronrod Rule */


/*           The abscissae and weights are given for the interval (-1,1) */
/*           because of symmetry only the positive abscisae and their */
/*           corresponding weights are given. */

/*           XGK    - abscissae of the 2N+1-point Kronrod rule: */
/*                    XGK(2), XGK(4), ...  N-point Gauss rule abscissae;*/ 
/*                    XGK(1), XGK(3), ...  abscissae optimally added */
/*                    to the N-point Gauss rule. */

/*           WGK    - weights of the 2N+1-point Kronrod rule. */
/*           WG     - weights of the N-point Gauss rule. */


/*     List of major variables */
/*           CENTER  - mid point of the interval */
/*           HFLGTH  - half-length of the interval */
/*           ABSCIS   - abscissae */
/*           RESLTG   - result of the N-point Gauss formula */
/*           RESLTK   - result of the 2N+1-point Kronrod formula */

    hflgth = (*b - *a) / 2;
    center = (*b + *a) / 2;

/*           compute the 2N+1-point Kronrod approximation to */
/*           the integral, and estimate the absolute error. */
    fc = trvfnd_(&center);
    resltg = fc * wg[0];
    resltk = fc * wgk[0];
    for (j = 1; j <= 11; ++j) {
	abscis = hflgth * xgk[j];
	d__1 = center - abscis;
	d__2 = center + abscis;
	funsum = trvfnd_(&d__1) + trvfnd_(&d__2);
	resltk += wgk[j] * funsum;
	if (j % 2 == 0) {
	    resltg += wg[j / 2] * funsum;
	}
    }
    ret_val = resltk * hflgth;
    *abserr = (d__1 = (resltk - resltg) * hflgth, abs(d__1)) * 3;
    return ret_val;
} /* krnrdd_ */


doublereal phid_(doublereal *z__)
{
    /* System generated locals */
    doublereal ret_val, d__1;
    /* Local variables */
    static doublereal zabs, p, expntl;

/*     Normal distribution probabilities accurate to 1.e-15. */
/*     Z = no. of standard deviations from the mean. */
/*     Based upon algorithm 5666 for the error function, from: */
/*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968 */
/*     Programmer: Alan Miller */
/*     Latest revision - 30 March 1986 */

    zabs = abs(*z__);

/*     |Z| > 37 */
    if (zabs > 37.) {
	p = 0.;
    } else {

/*     |Z| <= 37 */

/* Computing 2nd power */
	d__1 = zabs;
	expntl = exp(-(d__1 * d__1) / 2);

/*     |Z| < CUTOFF = 10/SQRT(2) */
	if (zabs < 7.071067811865475) {
	    p = expntl * ((((((zabs * .03526249659989109 + .7003830644436881) 
		    * zabs + 6.37396220353165) * zabs + 33.912866078383) * 
		    zabs + 112.0792914978709) * zabs + 221.2135961699311) * 
		    zabs + 220.2068679123761) / (((((((zabs * 
		    .08838834764831844 + 1.755667163182642) * zabs + 
		    16.06417757920695) * zabs + 86.78073220294608) * zabs + 
		    296.5642487796737) * zabs + 637.3336333788311) * zabs + 
		    793.8265125199484) * zabs + 440.4137358247522);

/*     |Z| >= CUTOFF. */
	} else {
	    p = expntl / (zabs + 1 / (zabs + 2 / (zabs + 3 / (zabs + 4 / (
		    zabs + .65))))) / 2.506628274631001;
	}
    }
    if (*z__ > 0.) {
	p = 1 - p;
    }
    ret_val = p;
    return ret_val;
} /* phid_ */


//Alternative gamma and tri and bivariate normal probabilities 

//Sub1
doublereal bvnorm_(doublereal *x, doublereal *y, doublereal *r__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;

/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*  Ensures that correct sign is used for parameters passed to Bivariate */
/*  normal integration routine. Alternatively, you can call BIVAR directly */
/*  using the opposite signs. Using BVNORM, you can write functions as you */
/*  see them in print... less confusing. */
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
    d__1 = -(*x);
    d__2 = -(*y);
    ret_val = bivar_(&d__1, &d__2, r__);
    return ret_val;
}

//Sub2
int bvn_(doublereal *x, doublereal *y, doublereal *r__, doublereal *p,
	doublereal *f, doublereal *px, doublereal *py, doublereal *p1, doublereal *f1, 
	doublereal *p2, doublereal *f2)
{
    /* System generated locals */
    doublereal d__1, d__2;
    /* Local variables */
    static doublereal g, h__;
    static doublereal dr, fax, fay;
    static doublereal pax, pay;

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/* usage: CALL BVN(z1,z2,rho,bvp,f,px,py,p1,f1,p2,f2) */
/* Returns several useful functions. */
/* P1 marginal cdf (argument 1) */
/* P2 marginal cdf (argument 2) */
/* P  bivariate normal cdf given R = rho */
/* PX */
/* PY */
/* P1 */
/* F1      various useful things etc....*/
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/ 
    dr = 1. / sqrt(1. - *r__ * *r__);
    *p1 = phi_(x, f1, &g, &h__);
    *p2 = phi_(y, f2, &g, &h__);
    d__1 = -(*x);
    d__2 = -(*y);
    *p = bivar_(&d__1, &d__2, r__);
    d__1 = (*x - *r__ * *y) * dr;
    pax = phi_(&d__1, &fax, &g, &h__);
    d__1 = (*y - *r__ * *x) * dr;
    pay = phi_(&d__1, &fay, &g, &h__);
    *px = *f1 * pay;
    *py = *f2 * pax;
    *f = *f1 * fay * dr;
    return 0;
}


//Sub 4 Common Block Declarations 
struct bvn01_1_ {
    doublereal w[15], x[15];};

#define bvn01_1 (*(struct bvn01_1_ *) &bvn01_)

struct dpreal_1_ {
    doublereal rl0, rl1;
    integer i0, i1;};

#define dpreal_1 (*(struct dpreal_1_ *) &dpreal_)

/* Initialized data */
static struct {
    doublereal e_1[30];
    } bvn01_ = { .05303709733976105, .11284582465517608, .15082452315872363, 
	    .16279133631194213, .15185641060466367, .12593625823209979, 
	    .094198393058496453, .0640788141334108, .0398456458245284, 
	    .02272413644209539, .01191223548930554, .005748310643806657, 
	    .00255593490701438, .001047812282114606, 3.96170048170894e-4, 
	    .02110687265306352, .11122304843701245, .27339875290117911, 
	    .50775546039766938, .8144213676108329, 1.193559990964792, 
	    1.645373297397144, 2.1701027938568, 2.7680303764366516, 
	    3.4394792198475525, 4.1848147744876557, 5.0044458955656313, 
	    5.8988261184898432, 6.8684550925062301, 7.9138801847749976 };

static struct {
    doublereal e_1[2];
    integer e_2[2];
    } dpreal_ = { 0., 1., 0, 1 };


doublereal bivar_(doublereal *x1, doublereal *x2, doublereal *r__)
{
    /* System generated locals */
    doublereal ret_val, d__1, d__2;
    /* Local variables */
    static doublereal a;
    static integer i__;
    static doublereal p, q, s, v, d1, d2, d3, d4, f1, f2, f3;

/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*Routine returns the approximation to the Bivariate normal. I takes two C*/
/*Z-scores and Rho (the correlation coefficient) as input and returns a  C*/
/*probability. Requires the PHI function below. Note: The arguments are  C*/
/*OPPOSITE of what is expected so, when this is used directly - switch   C*/
/* the signs of X1 and X2. */
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
    d1 = dpreal_1.rl0;
    if (*x1 < dpreal_1.rl0) {
	d1 = dpreal_1.rl1;
    }
    d2 = dpreal_1.rl0;
    if (*x2 < dpreal_1.rl0) {
	d2 = dpreal_1.rl1;
    }
    d3 = dpreal_1.rl1 - d1 * 2.;
    d4 = dpreal_1.rl1 - d2 * 2.;
    ret_val = dpreal_1.rl0;
    p = dpreal_1.rl1 / sqrt(dpreal_1.rl1 - *r__ * *r__);
    for (i__ = dpreal_1.i1; i__ <= 15; ++i__) {
	a = d4 * (*r__ * bvn01_1.x[i__ - 1] * d3 + *r__ * *x1 - *x2) * p;
/* L9: */
/* Computing 2nd power */
	d__1 = bvn01_1.x[i__ - 1] + d3 * *x1;
	ret_val += bvn01_1.w[i__ - 1] * 2. * phi_(&a, &f1, &f2, &f3) * exp(
		bvn01_1.x[i__ - 1] - d__1 * d__1 * .5) / 5.0132564549262001;
    }
    d__1 = -(*x2);
    d__2 = -(*x1);
    ret_val = ret_val * d3 * d4 - d1 * d2 + d1 * phi_(&d__1, &q, &v, &s) + d2 
	    * phi_(&d__2, &q, &v, &s);
    if (ret_val <= dpreal_1.rl0) {
	ret_val = 1e-8;
    }
    if (ret_val >= dpreal_1.rl1) {
	ret_val = .999999999;
    }
    return ret_val;
} 


//Sub 4 Common Block Declarations 
struct s4dpreal_1_ {
    doublereal rl0, rl1;
    integer i0, i1;};

#define s4dpreal_1 (*(struct s4dpreal_1_ *) &s4dpreal_)

struct n01_1_ {
    doublereal a1, a2, a3, a4, a5, u1, u2;};

#define n01_1 (*(struct n01_1_ *) &n01_)

/* Initialized data */
static struct {
    doublereal e_1[2];
    integer e_2[2];
    } s4dpreal_ = { 0., 1., 0, 1 };

static struct {
    doublereal e_1[7];
    } n01_ = { 1.330274429, 1.821255978, 1.781477937, .356563782, .31938153, 
	    .2316419, .39894228 };

/* Table of constant values */
static doublereal c_b2 = 7.5;

doublereal phi_(doublereal *az, doublereal *f, doublereal *g, doublereal *h__)
{
    /* System generated locals */
    doublereal ret_val, d__1;
    /* Local variables */
    static doublereal r__, t, z__;

/* *******************************************************************C*/
/*     ROUTINE RETURNS CUMULATIVE NORMAL OF Z IN PHI, DENSITY IN F,   C*/
/*     -F/(1.-PHI) IN G, F/PHI IN H                                   C*/
/* *******************************************************************C*/ 
    z__ = *az;
    if (abs(z__) > 7.5) {
	z__ = d_sign(&c_b2, &z__);
    }
/* Computing 2nd power */
    d__1 = z__;
    *f = n01_1.u2 * exp(-(d__1 * d__1) / 2.);
    t = s4dpreal_1.rl1 / (s4dpreal_1.rl1 + n01_1.u1 * abs(z__));
    r__ = ((((n01_1.a1 * t - n01_1.a2) * t + n01_1.a3) * t - n01_1.a4) * t + 
	    n01_1.a5) * t;
    ret_val = *f * r__;
    if (z__ > s4dpreal_1.rl0) {
	ret_val = s4dpreal_1.rl1 - ret_val;
    }
    *g = -(*f) / (s4dpreal_1.rl1 - ret_val);
    *h__ = *f / ret_val;
    return ret_val;
}


//Sub 5 Common Block Declarations 
struct gammfn_1_ {
    doublereal a[10], b[10], c__[10];};

#define gammfn_1 (*(struct gammfn_1_ *) &gammfn_)

struct s5dpreal_1_ {
    doublereal rl0, rl1;
    integer zero, one;};

#define s5dpreal_1 (*(struct s5dpreal_1_ *) &s5dpreal_)

/* Initialized data */
static struct {
    doublereal e_1[30];
    } gammfn_ = { .941785597796, .004415381325, .056850436816, -.004219835396,
	     .00132680818, -1.89302448e-4, 3.6069226e-5, -6.056604e-6, 
	    1.05491e-6, -1.75834e-7, -.019028540418, .491415393029, 
	    -.056815747821, .008357821226, -.001333232857, 2.20313282e-4, 
	    -3.7040211e-5, 6.283636e-6, -1.070343e-6, 1.77756e-7, 
	    1.035272238161, -.476782254386, .104882904204, -.022256271815, 
	    .004589049484, -9.24546035e-4, 1.82783463e-4, -3.5578802e-5, 
	    6.829026e-6, -1.254139e-6 };

static struct {
    doublereal e_1[2];
    integer e_2[2];
    } s5dpreal_ = { 0., 1., 0, 1 };


/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* Subroutine to calculate Gamma function and its first two derivatives C */
/* Accurate to 7 digits.  Returns for the value of P                    C */
/*   G0    =  Gamma(P)                                                  C */
/*   G1    =  dlogGamma(P)/dP        = Digamma                          C */
/*   G2    =  d**2logGamma(P)/dP**2  = Trigamma                         C */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
int xgamma_(doublereal *p, doublereal *g0, doublereal *g1, doublereal *g2, integer *ind)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal t, w, x, y, ydsum, ymult, ytsum;
    static integer na;
    static doublereal ap, bp, cp;
    static integer iw;
    static doublereal ap1, bp1, cp1, ap2, bp2, cp2, rl2;
    static integer int__, lvn;

/* ---------------------------------------------------------------- */
/* CONSTANTS FOR GAMMA FUNCTION ROUTINE */
/* ---------------------------------------------------------------- */
    rl2 = 2.;
    lvn = 11;
    *ind = s5dpreal_1.zero;
    if (abs(*p) < 1e-5 || *p >= 57.) {
	goto L5;
    }
    na = (integer) (*p) - s5dpreal_1.one;
    if (*p < s5dpreal_1.rl0) {
	na -= s5dpreal_1.one;
    }
    x = *p - (doublereal) ((real) na);
    int__ = s5dpreal_1.zero;
    if ((d__1 = x - s5dpreal_1.rl1, abs(d__1)) < 1e-9) {
	int__ = s5dpreal_1.one;
    }
    if (int__ == s5dpreal_1.one && *p < s5dpreal_1.rl0) {
	goto L5;
    }
    ymult = s5dpreal_1.rl1;
    ydsum = s5dpreal_1.rl0;
    ytsum = s5dpreal_1.rl0;
    y = *p;
    if (na == s5dpreal_1.zero) {
	goto L2;
    }
    iw = na / abs(na);
    w = (doublereal) ((real) iw);
    na = abs(na);
    if (iw == -s5dpreal_1.one) {
	y -= s5dpreal_1.rl1;
    }
    i__1 = na;
    for (i__ = s5dpreal_1.one; i__ <= i__1; ++i__) {
	y -= w;
	ymult *= pow_di(&y, &iw);
	ydsum += w / y;
/* L1: */
	ytsum += w / (y * y);
    }
L2:
    if (int__ == s5dpreal_1.one) {
	goto L4;
    }
    t = rl2 * x - 3.;
    ap1 = s5dpreal_1.rl0;
    ap = s5dpreal_1.rl0;
    bp1 = s5dpreal_1.rl0;
    bp = s5dpreal_1.rl0;
    cp1 = s5dpreal_1.rl0;
    cp = s5dpreal_1.rl0;
    for (i__ = s5dpreal_1.one; i__ <= 10; ++i__) {
	ap2 = ap1;
	ap1 = ap;
	ap = gammfn_1.a[lvn - i__ - 1] + rl2 * ap1 * t - ap2;
	bp2 = bp1;
	bp1 = bp;
	bp = gammfn_1.b[lvn - i__ - 1] + rl2 * bp1 * t - bp2;
	cp2 = cp1;
	cp1 = cp;
/* L3: */
	cp = gammfn_1.c__[lvn - i__ - 1] + rl2 * cp1 * t - cp2;
    }
    *g0 = (ap - ap1 * t) * ymult;
    *g1 = bp - bp1 * t + ydsum;
    *g2 = cp - cp1 * t - ytsum;
    return 0;
L4:
    *g0 = ymult;
    *g1 = ydsum - .5772156649;
    *g2 = 1.6449340668 - ytsum;
    return 0;
L5:
    *ind = s5dpreal_1.one;
    *g0 = s5dpreal_1.rl0;
    *g1 = s5dpreal_1.rl0;
    *g2 = s5dpreal_1.rl0;
    return 0;
} 


//Sub 6  Common Block Declarations 
struct s6gammfn_1_ {
    doublereal a[10], b[10], c__[10];};

#define s6gammfn_1 (*(struct s6gammfn_1_ *) &s6gammfn_)

struct s6dpreal_1_ {
    doublereal rl0, rl1;
    integer zero, one;};

#define s6dpreal_1 (*(struct s6dpreal_1_ *) &s6dpreal_)

/* Initialized data */
static struct {
    doublereal e_1[30];
    } s6gammfn_ = { .941785597796, .004415381325, .056850436816, -.004219835396,
	     .00132680818, -1.89302448e-4, 3.6069226e-5, -6.056604e-6, 
	    1.05491e-6, -1.75834e-7, -.019028540418, .491415393029, 
	    -.056815747821, .008357821226, -.001333232857, 2.20313282e-4, 
	    -3.7040211e-5, 6.283636e-6, -1.070343e-6, 1.77756e-7, 
	    1.035272238161, -.476782254386, .104882904204, -.022256271815, 
	    .004589049484, -9.24546035e-4, 1.82783463e-4, -3.5578802e-5, 
	    6.829026e-6, -1.254139e-6 };

static struct {
    doublereal e_1[2];
    integer e_2[2];
    } s6dpreal_ = { 0., 1., 0, 1 };


/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
/* FUNCTION   to calculate Gamma function                               C */
/* Accurate to 7 digits.  Returns for the value of P                    C */
/* Usage P = GAM(Z)                                                     C */
/* CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC */
doublereal gam_(doublereal *p)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Local variables */
    static integer i__;
    static doublereal t, w, x, y, ydsum, ymult, ytsum;
    static integer na;
    static doublereal ap;
    static integer iw;
    static doublereal ap1, ap2, rl2;
    static integer ind, int__, lvn;

/* ---------------------------------------------------------------- */
/* CONSTANTS FOR GAMMA FUNCTION ROUTINE */
/* ---------------------------------------------------------------- */
    rl2 = 2.;
    lvn = 11;
    ind = s6dpreal_1.zero;
    if (abs(*p) < 1e-5 || *p >= 57.) {
	goto L5;
    }
    na = (integer) (*p) - s6dpreal_1.one;
    if (*p < s6dpreal_1.rl0) {
	na -= s6dpreal_1.one;
    }
    x = *p - (doublereal) ((real) na);
    int__ = s6dpreal_1.zero;
    if ((d__1 = x - s6dpreal_1.rl1, abs(d__1)) < 1e-9) {
	int__ = s6dpreal_1.one;
    }
    if (int__ == s6dpreal_1.one && *p < s6dpreal_1.rl0) {
	goto L5;
    }
    ymult = s6dpreal_1.rl1;
    ydsum = s6dpreal_1.rl0;
    ytsum = s6dpreal_1.rl0;
    y = *p;
    if (na == s6dpreal_1.zero) {
	goto L2;
    }
    iw = na / abs(na);
    w = (doublereal) ((real) iw);
    na = abs(na);
    if (iw == -s6dpreal_1.one) {
	y -= s6dpreal_1.rl1;
    }
    i__1 = na;
    for (i__ = s6dpreal_1.one; i__ <= i__1; ++i__) {
	y -= w;
	ymult *= pow_di(&y, &iw);
	ydsum += w / y;
/* L1: */
	ytsum += w / (y * y);
    }
L2:
    if (int__ == s6dpreal_1.one) {
	goto L4;
    }
    t = rl2 * x - 3.;
    ap1 = s6dpreal_1.rl0;
    ap = s6dpreal_1.rl0;
    for (i__ = s6dpreal_1.one; i__ <= 10; ++i__) {
	ap2 = ap1;
	ap1 = ap;
	ap = s6gammfn_1.a[lvn - i__ - 1] + rl2 * ap1 * t - ap2;
/* L3: */
	ret_val = (ap - ap1 * t) * ymult;
    }
    return ret_val;
L4:
    ret_val = ymult;
    return ret_val;
L5:
    ind = s6dpreal_1.one;
    ret_val = s6dpreal_1.rl0;
    return ret_val;
}


int tnorm_(doublereal *z__, doublereal *x1, doublereal *x2, doublereal *x3, 
	doublereal *rho12, doublereal *rho31, doublereal *rho23)
{
    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal sxab[6], phix[6];
    static doublereal a[6], b[6], f, g, h__;
    static integer i__, j;
    static doublereal r__[6];
    static doublereal x[6], y[6], rr, xr[6], yr[3];
    static integer ind;
    static doublereal det;
    static doublereal txa[6], xiv[3], yrt[6];

/*     FOR THE TRIVARIATE NORMAL R.V. (X,Y,W) WITH ZERO MEANS, UNIT */
/*     VARIANCE AND CORR(X,Y)=RHO12,CORR(Y,W)RHO23,CORR(W,X)=RHO31,TNOR */
/*     USES THE METHOD OF SECTION 2 OF STECK (1958) TO COMPUTE Z= PR */
/*     (X<X1,Y<X2,W<X3) ACCURACY OF THE RESULT DEPENDS ON THE */
/*     SUBPROGRAMS T AND STK. FOR MORE DETAIL RE ACCURACY SEE DALEY */
/*     (1973) ANU STAT DEPT (IAS) REP. REF.: G.P. STECK (1958) ANN MATH */
/*     STATIS VOL 29,PP 780-800. */

/*     format is  call tnorm(p,1.97,-1.96,2.84,.1,.2,.3,999) */
/*         (  999 is alternate return label  ) */
/*        SUBROUTINE T AND STK USE CLOSED FORM APPROXIMATIONS */
/*        ACCURATE TO ABOUT .00025,  TIME IS ABOUT .003 */
    x[0] = *x1;
    x[1] = *x2;
    x[2] = *x3;
    for (i__ = 1; i__ <= 3; ++i__) {
	if ((d__1 = x[i__ - 1], abs(d__1)) >= (float)1e-15) {
	    goto L8;
	}
	x[i__ - 1] = (float)1e-15;
L8:
	;
    }
    r__[0] = *rho12;
    r__[1] = *rho23;
    r__[2] = *rho31;
    ind = 0;
/*     IND=0 FOR ALL X'X OF SAME SIGN, =1 OTHERWISE. IF X'S ARE OF */
/*     DIFFERENT SIGNS, THEN PERMUTE TO HAVE X(3) OF DIFF'T SIGN */
/*     AND CHANGE X(3),R(2), AND R(3). SIGNS. */
    if (x[0] != 0.) {
	goto L2;
    } else {
	goto L1;
    }
L1:
    if (x[1] * x[2] >= 0.) {
	goto L10;
    } else {
	goto L7;
    }
L2:
    if (x[0] * x[1] >= 0.) {
	goto L3;
    } else {
	goto L4;
    }
L3:
    if (x[0] * x[2] >= 0.) {
	goto L10;
    } else {
	goto L7;
    }
L4:
    if (x[0] * x[2] >= 0.) {
	goto L6;
    } else {
	goto L5;
    }
L5:
    swap_(x, &x[2]);
    swap_(r__, &r__[1]);
    goto L7;
L6:
    swap_(&x[1], &x[2]);
    swap_(r__, &r__[2]);
L7:
    ind = 1;
    x[2] = -x[2];
    r__[1] = -r__[1];
    r__[2] = -r__[2];
L10:
    for (i__ = 1; i__ <= 3; ++i__) {
/* L11: */
	y[i__ - 1] = (float)1. - r__[i__ - 1] * r__[i__ - 1];
    }
    det = y[0] + y[1] + y[2] - (float)2. + r__[0] * (float)2. * r__[1] * r__[
	    2];
    if (det > (float)1e-15) {
	goto L13;
    }
/*     IF(DET.LT.-1.E-4)GO TO 31 */
    if (det < (float)-1e-4) {
	return 1;
    }
    rr = (float).999999;
L13:
    det = (float)1. / sqrt(det);
    *z__ = (float)0.;
    yrt[0] = (float)1. / sqrt(y[0]);
    yrt[1] = (float)1. / sqrt(y[1]);
    yrt[2] = (float)1. / sqrt(y[2]);
    for (i__ = 1; i__ <= 3; ++i__) {
	j = i__ + 3;
	x[j - 1] = x[i__ - 1];
	xiv[i__ - 1] = (float)1. / x[i__ - 1];
	phix[i__ - 1] = phi_(&x[i__ - 1], &f, &g, &h__);
	phix[j - 1] = phix[i__ - 1];
	y[j - 1] = y[i__ - 1];
	yrt[j - 1] = yrt[i__ - 1];
	r__[j - 1] = r__[i__ - 1];
	yr[i__ - 1] = r__[i__ - 1] * r__[i__ + 1] - r__[i__];
	xr[i__ - 1] = x[i__] - x[i__ - 1] * r__[i__ - 1];
	if ((d__1 = xr[i__ - 1], abs(d__1)) > (float)1e-25) {
	    goto L14;
	}
	xr[i__ - 1] = (float)1e-25;
L14:
	xr[j - 1] = x[i__ + 1] - x[i__ - 1] * r__[i__ + 1];
	if ((d__1 = xr[j - 1], abs(d__1)) >= (float)1e-25) {
	    goto L15;
	}
	xr[j - 1] = (float)1e-25;
L15:
	a[i__ - 1] = xr[i__ - 1] * xiv[i__ - 1] * yrt[i__ - 1];
	a[j - 1] = xr[j - 1] * xiv[i__ - 1] * yrt[i__ + 1];
	b[i__ - 1] = det * (yr[i__ - 1] + y[i__ - 1] * xr[j - 1] / xr[i__ - 1]
		);
	b[j - 1] = det * (yr[i__ - 1] + y[i__ + 1] * xr[i__ - 1] / xr[j - 1]);
	*z__ += ((float)1. - del_(&a[i__ - 1], &a[j - 1])) * phix[i__ - 1];
/* L16: */
    }
    *z__ *= (float).5;
    for (i__ = 1; i__ <= 6; ++i__) {
	txa[i__ - 1] = t_(&x[i__ - 1], &a[i__ - 1], &phix[i__ - 1]);
	sxab[i__ - 1] = s_(&x[i__ - 1], &a[i__ - 1], &b[i__ - 1], &phix[i__ - 
		1], &txa[i__ - 1]);
	*z__ = *z__ - txa[i__ - 1] * (float).5 - sxab[i__ - 1];
/* L17: */
    }
    if (ind == 0) {
	return 0;
    }
    *z__ = (phix[0] + phix[4] - del_(x, &x[1])) * (float).5 - txa[0] - txa[4] 
	    - *z__;
    return 0;
} /* tnorm_ */


int swap_(doublereal *x, doublereal *y)
{
    static doublereal z__;
    z__ = *x;
    *x = *y;
    *y = z__;
    return 0;
} 

doublereal del_(doublereal *x, doublereal *y)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    if ((d__1 = *x * *y) < 0.) {
	goto L1;
    } else if (d__1 == 0) {
	goto L3;
    } else {
	goto L2;
    }
L1:
    ret_val = (float)1.;
    return ret_val;
L2:
    ret_val = (float)0.;
    return ret_val;
L3:
    if (*x + *y >= 0.) {
	goto L2;
    } else {
	goto L1;
    }
} 

doublereal s_(doublereal *x, doublereal *a, doublereal *b, doublereal *phix, doublereal *txa)
{
	/*COMPUTES STECK'S (1958 ANN MATH STAT) S-FN WITH MAX ERROR E-9, */
	/*THOUGH COMPUTER ROUNDING ERRORS MAY REDUCE THIS ACCURACY. CALLS */
	/*FUNCTION SUBPROGRAMS STK(X,A,B) IN WHICH DABS(B) < 1.1, PHI(X), AND */
	/*T(X,A,PHIX) NOTE THAT PHIX = PHI(X) AND TXA = T(X,A,PHIX). */

    /* System generated locals */
    doublereal ret_val, d__1;
    /* Local variables */
    static doublereal f, g, h__;
    static doublereal aa, ba, ab, ta, xa, xf, phixab, xab, xfa, biv;
    ba = abs(*b);
    if (ba > (float)1.1) {
	goto L1;
    }
    ret_val = stk_(x, a, b);
    return ret_val;
L1:
    ta = abs(*txa);
    aa = abs(*a);
    biv = (float)1. / ba;
    xf = (d__1 = *phix - (float).5, abs(d__1));
    xa = abs(*x);
    xab = xa * aa * ba;
    phixab = phi_(&xab, &f, &g, &h__);
    xfa = phixab - (float).5;
    if (aa <= (float)1.) {
	goto L2;
    }
    d__1 = (float)1. / aa;
    ret_val = xfa * ta - xf * t_(&xab, &biv, &phixab) + (xf - xfa) * (float)
	    .25 + stk_(&xab, &biv, &d__1);
    goto L3;
L2:
    ab = aa * ba;
    d__1 = (float)1. / ab;
    ret_val = (xf + (float).5) * (float).25 + xfa * ta - stk_(&xab, &d__1, &
	    aa) - stk_(&xa, &ab, &biv);
L3:
    if (*x >= (float)0.) {
	goto L4;
    }
    ret_val = atan(ba / sqrt(*a * *a * (*b * *b + (float)1.) + (float)1.)) * (
	    float).1591549 - ret_val;
L4:
    if (*b >= (float)0.) {
	return ret_val;
    }
    ret_val = -ret_val;
    return ret_val;
} 


doublereal t_(doublereal *z__, doublereal *a, doublereal *phix)
{
	/* OWENS T DOUBLE PRECISION FUNCTION ACCURATE TO .00005 */
	/* PHIX IS A DUMMY ARGUMENT */

    /* System generated locals */
    doublereal ret_val, d__1;
    /* Local variables */
    static doublereal f, g, h__, theta, a1, f1, x1, ff1, gg1, gg2, x2a;
    x1 = *z__;
    a1 = *a;
    ff1 = (float)0.;
    f1 = (float)1.;
    ret_val = (float)0.;
    if (*a >= (float)0.) {
	goto L10;
    }
    a1 = -a1;
    f1 = (float)-1.;
L10:
    if (a1 <= (float)1.) {
	goto L15;
    }
    x1 = abs(x1);
    gg1 = phi_(&x1, &f, &g, &h__);
    x1 = a1 * x1;
    a1 = (float)1. / a1;
    gg2 = phi_(&x1, &f, &g, &h__);
    ff1 = (gg1 + gg2) * (float).5 - gg1 * gg2;
    ff1 *= f1;
    f1 = -f1;
L15:
    theta = atan(a1);
    x2a = x1 * x1;
    if (x2a > (float)150.) {
	goto L20;
    }
    ret_val = f1 * theta * (float).15915494 * exp(x2a * (float)-.5 * a1 / theta);
	/* Computing 4th power */
    d__1 = x1 * a1, d__1 *= d__1;
    ret_val *= d__1 * d__1 * (float).00868 + (float)1.;
    ret_val += ff1;
L20:
    return ret_val;
} 


doublereal stk_(doublereal *x, doublereal *a, doublereal *b)
{
    /* System generated locals */
    doublereal ret_val, d__1;
    /* Local variables */
    static doublereal f, g, h__, z__, z1, z2;

	/*ACCURATE TO .0001 */
    z__ = *a * *a + (float)1.;
	/* Computing 2nd power */
    d__1 = *a * (float).5 * *b;
    z1 = sqrt(z__ + d__1 * d__1);
	/* Computing 2nd power */
    d__1 = *a * *b;
    z2 = *b / sqrt(z__ + d__1 * d__1);
    d__1 = *x * z1;
    ret_val = phi_(&d__1, &f, &g, &h__) * (float).1591549 * atan(z2);
    return ret_val;
} 


doublereal den2_(doublereal *z1, doublereal *z2, doublereal *r__)
{
    /* Initialized data */
    static doublereal twopi = 6.28318530718;
    /* System generated locals */
    doublereal ret_val;
    /* Local variables */
    static doublereal term1, term2, onemr2;

	/* FUNCTION TO COMPUTE STANDARD BINORMAL DENSITY */
    onemr2 = 1. - *r__ * *r__;
    term1 = 1. / (twopi * sqrt(onemr2));
    term2 = -(*z1 * *z1 - *r__ * 2. * *z1 * *z2 + *z2 * *z2) / (onemr2 * 2.);
    ret_val = term1 * exp(term2);
    return ret_val;
} 


/* Main program
void main()
{
//Get input parameters
double lim[3], sig[3];
int opt;
float x1, x2, x3, sig1, sig2, sig3, z1;
doublereal xx1, xx2, xx3, ssig1, ssig2, ssig3, zz1;

	opt = -1;

	while (opt != 0) {
// Get option 
	printf ("1- Approximation  2 Integration  - 0 Exit \n");
	scanf("%d",&opt);
	if (opt != 0){
	printf ("\nX1 ");
	scanf("%f",&x1);
	printf ("\nX2 ");
	scanf("%f",&x2);
	printf ("\nX3 ");
	scanf("%f",&x3);
	printf ("\nSig1 ");
	scanf("%f",&sig1);
	printf ("\nSig2 ");
	scanf("%f",&sig2);
	printf ("\nSig3 ");
	scanf("%f",&sig3);
	printf ("\n");
	}

xx1=x1; xx2=x2; xx3=x3;
ssig1=sig1; ssig2=sig2; ssig3=sig3;

if (opt == 1) {
	zz1=0.0;
	tnorm_(&zz1, &xx1, &xx2, &xx3, &ssig1, &ssig2, &ssig3);}
if (opt == 2) {
	zz1=0.0;
	xx1=x1; xx2=x2; xx3=x3;
	lim[0]=xx1; lim[1]=xx2; lim[2]=xx3;
	sig[0]=ssig1; sig[1]=ssig2; sig[2]=ssig3;
    zz1=tvnd_(lim, sig);}
if (opt == 0)
    return;

    z1 = (float)zz1;
    printf ("\nProbability is:   %e \n", z1);
	printf ("\n\n");

	}  //While

} 
 */



