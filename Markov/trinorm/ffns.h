typedef double doublereal;
typedef long int integer;
typedef float real;
typedef long int logical;

#define TRUE_ (1)
#define FALSE_ (0)

#define log10e 0.43429448190325182765
#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)

double d_nint(doublereal *x);
double d_int(doublereal *x);
double pow_di(doublereal *ap, integer *bp);
double d_sign(doublereal *a, doublereal *b);
double pow_dd(doublereal *ap, doublereal *bp);
double d_lg10(doublereal *x);


doublereal bvn_(real *sh, real *sk, real *r__);
doublereal phi_(real *x);
doublereal bvnd_(doublereal *dh, doublereal *dk, doublereal *r__);
doublereal tvnd_(doublereal *limit, doublereal *sigma);
doublereal trvfnd_(doublereal *t);
doublereal adoned_(doublereal *a, doublereal *b, doublereal *tol);
doublereal krnrdd_(doublereal *a, doublereal *b, doublereal *abserr);
doublereal phid_(doublereal *z__);


doublereal bvnorm_(doublereal *x, doublereal *y, doublereal *r__);
int bvn_(doublereal *x, doublereal *y, doublereal *r__, doublereal *p,
	doublereal *f, doublereal *px, doublereal *py, doublereal *p1, doublereal *f1, 
	doublereal *p2, doublereal *f2);
doublereal bivar_(doublereal *x1, doublereal *x2, doublereal *r__);
doublereal phi_(doublereal *az, doublereal *f, doublereal *g, doublereal *h__);
int xgamma_(doublereal *p, doublereal *g0, doublereal *g1, doublereal *g2, integer *ind);
doublereal gam_(doublereal *p);
int tnorm_(doublereal *z__, doublereal *x1, doublereal *x2, doublereal *x3, 
	doublereal *rho12, doublereal *rho31, doublereal *rho23);
int swap_(doublereal *x, doublereal *y);
doublereal del_(doublereal *x, doublereal *y);
doublereal s_(doublereal *x, doublereal *a, doublereal *b, doublereal *phix, doublereal *txa);
doublereal t_(doublereal *z__, doublereal *a, doublereal *phix);
doublereal stk_(doublereal *x, doublereal *a, doublereal *b);
doublereal den2_(doublereal *z1, doublereal *z2, doublereal *r__);
