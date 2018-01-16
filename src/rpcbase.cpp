#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "rpcbase.h"

#include <gdal_priv.h>
#include <gdal_alg.h>
#include "PrjLog.hpp"

#ifndef WIN32

#ifndef _strlwr
#define _strlwr
inline char* strlwr(char *str){
    char *orig = str; char d = 'a' - 'A';
    for (; *str != '\0'; str++){ if (*str >= 'A' && *str <= 'Z') *str = *str + d; }
    return orig;
}
#endif

inline char* strupr(char *str){
    char *orig = str; char d = 'a' - 'A';
    for (; *str != '\0'; str++){ if (*str >= 'a' && *str <= 'z') *str = *str - d; }
    return orig;
}

inline int stricmp(const char *s1,const char *s2){
    char strT1[512], strT2[512];
    strcpy(strT1, s1); strcpy(strT2, s2);
    strlwr(strT1); strlwr(strT2);
    return strcmp(strT1, strT2);
}

inline int strnicmp(const char *s1,const char *s2, int n){
    char strT1[512], strT2[512];
    strcpy(strT1, s1); strcpy(strT2, s2);
    strlwr(strT1); strlwr(strT2);
    return strncmp(strT1, strT2, n);
}
#endif // WIN32

static void _ati_ldltban1(double *a, double *d, double *l, int n, int wide){
    int i, j, k, kk, km, m;
    double *ao, *aa, *co, *c;

    m = wide*(2 * n + 1 - wide) / 2;
    c = new double[(m - wide)];

    ao = a; co = c; a += wide;
    for (i = 0; i < m - wide; i++){
        *c++ = *a++;
    }
    c = co; a = ao;

    kk = 0;

    for (k = 1; k < n; k++) {
        if (k < n - wide + 2) kk = wide - 1;
        else kk--;

        *d = *a++; aa = a;  a += kk;

        if (k < n - wide + 1) km = wide;
        else km = n - k + 1;

        for (i = 1; i < kk + 1; i++) {
            *l = *aa++ / *d;
            for (j = 0; j<kk - i + 1; j++) *(a + j) -= *l * *(aa + j - 1);
            l++;

            if (k + i>n - wide + 1) km--;
            a += km;
        }

        a = aa; d++;
        if (k == n - 1)  *d = *a;
    }

    a = ao;  a += wide;
    for (i = 0; i < m - wide; i++){
        *a++ = *c++;
    }
    c = co; delete[]c;
};

static void _ati_ldltban2(double *l, double *d, double *b, double *x, int n, int wide){
    int i, j, kk, m;
    double *bo, *lo, *xx;
    double *bb, *bbo;
    kk = 0;

    bb = new double[n];
    bbo = bb;

    bo = b; lo = l;

    for (i = 0; i < n; i++)*bb++ = *b++;
    b = bo;  bb = bbo;
    m = wide*(2 * n + 1 - wide) / 2;

    for (i = 1; i < n; i++){
        if (i < n - wide + 2) kk = wide;
        else kk--;

        b = bo + i;
        for (j = 1; j < kk; j++){
            *b++ -= *(b - j) * *l++;
        }
    }

    kk = 0;
    b = bo + n - 1;  l = lo + m - n - 1;
    x += n - 1;  xx = x;  d += n - 1;

    *x-- = (*b--) / (*d--);

    for (i = 1; i < n; i++)  {
        if (i < wide) kk++;
        else { kk = wide - 1;  xx--; }

        *x = *b-- / *d--;
        for (j = 1; j < kk + 1; j++){
            *x -= *l-- * *(xx - j + 1);
        }
        x--;
    }

    b = bo;
    for (i = 0; i < n; i++){
        *b++ = *bb++;
    }
    bb = bbo; delete[]bb;
};

// construct the normalize equations
static void _ati_nrml(double *aa, int n, double bb, double *a, double *b){
    register int  i, j;
    double *a0 = a;

    for (i = 0; i < n; i++){
        for (j = 0; j < n - i; j++){
            *a += *aa * *(aa + j);
            a++;
        }
        *b += *aa * bb;
        b++; aa++;
    }

    a = a0;
};

// construct the normalize equations
static void _ati_pnrml(double *aa, int n, double bb, double *a, double *b, double p){
    register int  i, j;
    double *a0 = a;

    for (i = 0; i < n; i++){
        for (j = 0; j < n - i; j++){
            *a += ((*aa)*(*(aa + j))*p);
            a++;
        }
        *b += ((*aa)*bb*p);
        b++; aa++;
    }

    a = a0;
};

// solve the normalized equation
static void _ati_solve(double*a, double*b, double*x, int n, int wide){
    int    m;
    double *d, *l;

    m = n*(n + 1) / 2;
    d = new double[n];
    l = new double[m - n];

    memset(d, 0, n*sizeof(double));
    memset(l, 0, (m - n)*sizeof(double));

    _ati_ldltban1(a, d, l, n, wide);
    _ati_ldltban2(l, d, b, x, n, wide);

    delete[]d;
    delete[]l;
};

// calculate the polynominal coeffcients of P, L and H
void _ati_rpc_plh_coef(double P, double L, double H, double* coef){
    coef[0] = 1.0;
    coef[1] = L;
    coef[2] = P;
    coef[3] = H;
    coef[4] = L*P;
    coef[5] = L*H;
    coef[6] = P*H;
    coef[7] = L*L;
    coef[8] = P*P;
    coef[9] = H*H;
    coef[10] = P*L*H;
    coef[11] = L*L*L;
    coef[12] = L*P*P;
    coef[13] = L*H*H;
    coef[14] = L*L*P;
    coef[15] = P*P*P;
    coef[16] = P*H*H;
    coef[17] = L*L*H;
    coef[18] = P*P*H;
    coef[19] = H*H*H;
};

// calculating 20 coefficients of dF/dP
void _ati_differential_p(double P, double L, double H, double* coef){
    memset(coef, 0, 20 * sizeof(double));
    coef[2] = 1;
    coef[4] = L;
    coef[6] = H;
    coef[8] = 2 * P;
    coef[10] = L*H;
    coef[12] = 2 * L*P;
    coef[14] = L*L;
    coef[15] = 3 * P*P;
    coef[16] = H*H;
    coef[18] = 2 * P*H;
};

// calculating 20 coefficients of dF/dL
void _ati_differential_l(double P, double L, double H, double* coef){
    memset(coef, 0, 20 * sizeof(double));
    coef[1] = 1;
    coef[4] = P;
    coef[5] = H;
    coef[7] = 2 * L;
    coef[10] = P*H;
    coef[11] = 3 * L*L;
    coef[12] = P*P;
    coef[13] = H*H;
    coef[14] = 2 * L*P;
    coef[17] = 2 * L*H;
};

// calculating 20 coefficients of dF/dH
void _ati_differential_h(double P, double L, double H, double* coef){
    memset(coef, 0, 20 * sizeof(double));
    coef[3] = 1;
    coef[5] = L;
    coef[6] = P;
    coef[9] = 2 * H;
    coef[10] = P*L;
    coef[13] = 2 * L*H;
    coef[16] = 2 * P*H;
    coef[17] = L*L;
    coef[18] = P*P;
    coef[19] = 3 * H*H;
};

// zoom the vector by mutiple
static void _ati_rpc_numeric_product(double* vector, int num, double multiple){
    for (int i = 0; i < num; i++){
        vector[i] *= multiple;
    }
};

// calculate the scalar product of two vectors
static double _ati_rpc_scalar_product(double* vector1, double* vector2, int number){
    double product = 0.0;
    for (int i = 0; i < number; i++)
        product += (vector1[i] * vector2[i]);
    return product;
};

// project the (lat, lon, alt) to image through rpc parameters
// the reference of output image coordinates is on the left-upper corner
static void _rpc_geoedtic_to_pxy(RpcPara* pRpc, double lat, double lon, double alt, double* x, double* y){
    // normalize the lat, lon, and alt
    double P = (lat - pRpc->lat_off) / (pRpc->lat_scale);
    double L = (lon - pRpc->long_off) / (pRpc->long_scale);
    double H = (alt - pRpc->height_off) / (pRpc->height_scale);

    // calculate the PLH coefficients
    double coef[20]; _ati_rpc_plh_coef(P, L, H, coef);

    // calculate the normalizing photo coordinates
    double Y = _ati_rpc_scalar_product(pRpc->c, coef, 20) / _ati_rpc_scalar_product(pRpc->d, coef, 20);
    double X = _ati_rpc_scalar_product(pRpc->a, coef, 20) / _ati_rpc_scalar_product(pRpc->b, coef, 20);

    // export the photo coordinates (the reference point is on the left-upper corner)
    (*y) = Y*pRpc->line_scale + pRpc->line_off;
    (*x) = X*pRpc->samp_scale + pRpc->samp_off;
};

// Forward Intersection
static void _rpc_pxy_to_geoedtic(RpcPara* lrpc, double lx, double ly,
    RpcPara* rrpc, double rx, double ry,
    double *lat, double *lon, double *alt)
{
    // assign the initial value for output
    *lat = *lon = *alt = 0.0;

    // normalize the photo coordinates
    RpcPara* p = lrpc;
    lx = (lx - p->samp_off) / (p->samp_scale);
    ly = (ly - p->line_off) / (p->line_scale);

    p = rrpc;
    rx = (rx - p->samp_off) / (p->samp_scale);
    ry = (ry - p->line_off) / (p->line_scale);

    // use the constant and linear coefficients of a, b, c, d
    // to calculate the initial geoedtic coordinates
    double latitude, longitude, height, F0, G0;
    latitude = lrpc->lat_off;
    longitude = lrpc->long_off;
    height = lrpc->height_off;

    // iterative calculat the geoedtic coordinates
    double x, y, P, L, H, A[3], B, ATA[9], ATB[3], X[3];
    double v[20], dp[20], dl[20], dh[20], f[20];

    int k, iters = 0;
    do{
        // normalize the error equations
        memset(ATA, 0, 9 * sizeof(double));
        memset(ATB, 0, 3 * sizeof(double));

        /////////////////////////////////////////////
        // Left Image
        p = lrpc; x = lx; y = ly;
        // calculate the normalize geoedtic coordinates
        P = (latitude - p->lat_off) / (p->lat_scale);
        L = (longitude - p->long_off) / (p->long_scale);
        H = (height - p->height_off) / (p->height_scale);
        // calculate the polynominal coeffcients of P, L and H
        _ati_rpc_plh_coef(P, L, H, f);
        // calculate differential coeffcients of P, L and H
        _ati_differential_p(P, L, H, dp);
        _ati_differential_l(P, L, H, dl);
        _ati_differential_h(P, L, H, dh);

        // normalize the error equation of x coordinate
        F0 = _ati_rpc_scalar_product(p->a, f, 20);
        G0 = _ati_rpc_scalar_product(p->b, f, 20);
        for (k = 0; k < 20; k++)
            v[k] = G0*(p->a[k]) - F0*(p->b[k]);

        A[0] = _ati_rpc_scalar_product(v, dp, 20) / (p->lat_scale) / G0 / G0;
        A[1] = _ati_rpc_scalar_product(v, dl, 20) / (p->long_scale) / G0 / G0;
        A[2] = _ati_rpc_scalar_product(v, dh, 20) / (p->height_scale) / G0 / G0;
        B = x - F0 / G0;
        _ati_pnrml(A, 3, B, ATA, ATB, 1.0);

        // normalize the error equation of y coordinates
        F0 = _ati_rpc_scalar_product(p->c, f, 20);
        G0 = _ati_rpc_scalar_product(p->d, f, 20);
        for (k = 0; k < 20; k++)
            v[k] = G0*(p->c[k]) - F0*(p->d[k]);

        A[0] = _ati_rpc_scalar_product(v, dp, 20) / (p->lat_scale) / G0 / G0;
        A[1] = _ati_rpc_scalar_product(v, dl, 20) / (p->long_scale) / G0 / G0;
        A[2] = _ati_rpc_scalar_product(v, dh, 20) / (p->height_scale) / G0 / G0;
        B = y - F0 / G0;
        _ati_pnrml(A, 3, B, ATA, ATB, 1.0);

        /////////////////////////////////////////////
        // Right Image
        p = rrpc; x = rx; y = ry;
        // calculate the normalize geoedtic coordinates
        P = (latitude - p->lat_off) / (p->lat_scale);
        L = (longitude - p->long_off) / (p->long_scale);
        H = (height - p->height_off) / (p->height_scale);
        // calculate the polynominal coeffcients of P, L and H
        _ati_rpc_plh_coef(P, L, H, f);
        // calculate differential coeffcients of P, L and H
        _ati_differential_p(P, L, H, dp);
        _ati_differential_l(P, L, H, dl);
        _ati_differential_h(P, L, H, dh);

        // normalize the error equation of x coordinate
        F0 = _ati_rpc_scalar_product(p->a, f, 20);
        G0 = _ati_rpc_scalar_product(p->b, f, 20);
        for (k = 0; k < 20; k++)
            v[k] = G0*(p->a[k]) - F0*(p->b[k]);

        A[0] = _ati_rpc_scalar_product(v, dp, 20) / (p->lat_scale) / G0 / G0;
        A[1] = _ati_rpc_scalar_product(v, dl, 20) / (p->long_scale) / G0 / G0;
        A[2] = _ati_rpc_scalar_product(v, dh, 20) / (p->height_scale) / G0 / G0;
        B = x - F0 / G0;
        _ati_pnrml(A, 3, B, ATA, ATB, 1.0);

        // normalize the error equation of y coordinates
        F0 = _ati_rpc_scalar_product(p->c, f, 20);
        G0 = _ati_rpc_scalar_product(p->d, f, 20);
        for (k = 0; k < 20; k++)
            v[k] = G0*(p->c[k]) - F0*(p->d[k]);

        A[0] = _ati_rpc_scalar_product(v, dp, 20) / (p->lat_scale) / G0 / G0;
        A[1] = _ati_rpc_scalar_product(v, dl, 20) / (p->long_scale) / G0 / G0;
        A[2] = _ati_rpc_scalar_product(v, dh, 20) / (p->height_scale) / G0 / G0;
        B = y - F0 / G0;
        _ati_pnrml(A, 3, B, ATA, ATB, 1.0);

        // solve the normal equation to get the variance of unknowns
        memset(X, 0, 3 * sizeof(double));
        _ati_solve(ATA, ATB, X, 3, 3);
        latitude += X[0];
        longitude += X[1];
        height += X[2];

        // if iterative number reaches 50, break the iteration
        iters++;
        if (iters >= 25) break;
    } while (fabs(X[0]) > 1.0e-12 || fabs(X[1]) > 1.0e-12 || fabs(X[2]) > 1.0e-10);

    // output the geoedtic coordinates
    *lat = latitude;
    *lon = longitude;
    *alt = height;
}

static double _rpc_pxy_Z_to_geoedtic(RpcPara* pRpc, double sx, double sy, double alt0,
    double *lat, double *lon, double *alt)
{
    double lon0 = pRpc->long_off;
    double lat0 = pRpc->lat_off; int iter = 0;
    double sx0, sy0, sx1, sy1, sx2, sy2, dis, dx, dy;
    do{
        _rpc_geoedtic_to_pxy(pRpc, lat0, lon0, alt0, &sx0, &sy0);
        dis = sqrt((sx0 - sx)*(sx0 - sx) + (sy0 - sy)*(sy0 - sy));
        if (dis < 0.005) break;

        _rpc_geoedtic_to_pxy(pRpc, lat0, lon0 + 1.0e-06, alt0, &sx1, &sy1);
        _rpc_geoedtic_to_pxy(pRpc, lat0 + 1.0e-06, lon0, alt0, &sx2, &sy2);
        dx = ((sx - sx0)*(sx1 - sx0) + (sy - sy0)*(sy1 - sy0)) / ((sx1 - sx0)*(sx1 - sx0) + (sy1 - sy0)*(sy1 - sy0));
        dy = ((sx - sx0)*(sx2 - sx0) + (sy - sy0)*(sy2 - sy0)) / ((sx2 - sx0)*(sx2 - sx0) + (sy2 - sy0)*(sy2 - sy0));
        lon0 += (dx*1.0e-06);
        lat0 += (dy*1.0e-06);
        iter++;
    } while (iter <= 20);
    *lon = lon0; *lat = lat0; *alt = alt0;
    return dis;
}


// forward intersect to calculate the lat, lon and alt;
// this is a reductive version with no gross error detection,
// iterative power and post variance
static void _rpc_pxy_to_geoedtic(RpcPara** pRpc, double *pIx, double *pIy, int ptSum,
    double *lat, double *lon, double *alt)
{
    int i, k, n;
    // assign the initial value for output
    *lat = *lon = *alt = 0.0;
    // at least two observations
    if (ptSum > 2)
    {
        // normalize the photo coordinates
        RpcPara* p;
        for (i = 0; i < ptSum; i++){
            p = pRpc[i];
            pIx[i] = (pIx[i] - p->samp_off) / (p->samp_scale);
            pIy[i] = (pIy[i] - p->line_off) / (p->line_scale);
        }

        // use the constant and linear coefficients of a, b, c, d
        // to calculate the initial geoedtic coordinates
        double latitude, longitude, height, F0, G0;
        latitude = pRpc[0]->lat_off;
        longitude = pRpc[0]->long_off;
        height = pRpc[0]->height_off;

        // iterative calculat the geoedtic coordinates
        double x, y, P, L, H, A[3], B, ATA[9], ATB[3], X[3];
        double v[20], dp[20], dl[20], dh[20], f[20];

        int iters = 0;
        do{
            // normalize the error equations
            memset(ATA, 0, 9 * sizeof(double));
            memset(ATB, 0, 3 * sizeof(double));
            for (i = 0, n = 0; i < ptSum; i++)
            {
                p = pRpc[i];
                x = pIx[i];
                y = pIy[i];
                // calculate the normalize geoedtic coordinates
                P = (latitude - p->lat_off) / (p->lat_scale);
                L = (longitude - p->long_off) / (p->long_scale);
                H = (height - p->height_off) / (p->height_scale);
                // calculate the polynominal coeffcients of P, L and H
                _ati_rpc_plh_coef(P, L, H, f);
                // calculate differential coeffcients of P, L and H
                _ati_differential_p(P, L, H, dp);
                _ati_differential_l(P, L, H, dl);
                _ati_differential_h(P, L, H, dh);

                // normalize the error equation of x coordinate
                F0 = _ati_rpc_scalar_product(p->a, f, 20);
                G0 = _ati_rpc_scalar_product(p->b, f, 20);
                for (k = 0; k < 20; k++)
                    v[k] = G0*(p->a[k]) - F0*(p->b[k]);

                A[0] = _ati_rpc_scalar_product(v, dp, 20) / (p->lat_scale) / G0 / G0;
                A[1] = _ati_rpc_scalar_product(v, dl, 20) / (p->long_scale) / G0 / G0;
                A[2] = _ati_rpc_scalar_product(v, dh, 20) / (p->height_scale) / G0 / G0;
                B = x - F0 / G0;
                _ati_pnrml(A, 3, B, ATA, ATB, 1.0);

                // normalize the error equation of y coordinates
                F0 = _ati_rpc_scalar_product(p->c, f, 20);
                G0 = _ati_rpc_scalar_product(p->d, f, 20);
                for (k = 0; k < 20; k++)
                    v[k] = G0*(p->c[k]) - F0*(p->d[k]);

                A[0] = _ati_rpc_scalar_product(v, dp, 20) / (p->lat_scale) / G0 / G0;
                A[1] = _ati_rpc_scalar_product(v, dl, 20) / (p->long_scale) / G0 / G0;
                A[2] = _ati_rpc_scalar_product(v, dh, 20) / (p->height_scale) / G0 / G0;
                B = y - F0 / G0;
                _ati_pnrml(A, 3, B, ATA, ATB, 1.0);

                // accumulate the number of error quations
                n += 2;
            }
            if (n < 3) return;

            // solve the normal equation to get the variance of unknowns
            memset(X, 0, 3 * sizeof(double));
            _ati_solve(ATA, ATB, X, 3, 3);
            latitude += X[0];
            longitude += X[1];
            height += X[2];

            // if iterative number reaches 50, break the iteration
            iters++;
            if (iters >= 50) break;
        } while (fabs(X[0]) > 1.0e-12 || fabs(X[1]) > 1.0e-12 || fabs(X[2]) > 1.0e-10);

        // output the geoedtic coordinates
        *lat = latitude;
        *lon = longitude;
        *alt = height;
    }
}

static void _ati_rpc_off_scale(const double* p, int num, double& OFF, double& SCALE){
    double max, min;
    max = min = OFF = p[0];
    for (int i = 1; i < num; i++){
        if (max < p[i]) max = p[i];
        if (min > p[i]) min = p[i];
        OFF += p[i];
    }
    OFF /= num;

    double MAX = fabs(max - OFF);
    double MIN = fabs(min - OFF);
    SCALE = max - min;
};

// calculate 59 RPC parameters through the input 2D points and 3D points
static void _rpc_calc_rpc_59(double* samp, double* line, double* lat, double* lon, double* alt, int num, RpcPara* rpcPara, double& sx, double& sy)
{
    memset(rpcPara, 0, sizeof(RpcPara));
    _ati_rpc_off_scale(samp, num, rpcPara->samp_off, rpcPara->samp_scale);
    _ati_rpc_off_scale(line, num, rpcPara->line_off, rpcPara->line_scale);
    _ati_rpc_off_scale(lat, num, rpcPara->lat_off, rpcPara->lat_scale);
    _ati_rpc_off_scale(lon, num, rpcPara->long_off, rpcPara->long_scale);
    _ati_rpc_off_scale(alt, num, rpcPara->height_off, rpcPara->height_scale);
    rpcPara->samp_scale *= 2.0;
    rpcPara->line_scale *= 2.0;
    rpcPara->lat_scale *= 2.0;
    rpcPara->long_scale *= 2.0;

    double a[59], ata[59 * 59], atl[59], x[59];
    double P, L, H, X, Y, A[20], B[20], p, h, maxv;

    int i, n, k, iter = 0;
    do{
        memset(ata, 0, 59 * 59 * sizeof(double));
        memset(atl, 0, 59 * sizeof(double));
        memset(a, 0, 59 * sizeof(double));
        memset(x, 0, 59 * sizeof(double));

        for (i = 0; i < num; i++){
            X = (samp[i] - rpcPara->samp_off) / rpcPara->samp_scale;
            Y = (line[i] - rpcPara->line_off) / rpcPara->line_scale;
            P = (lat[i] - rpcPara->lat_off) / rpcPara->lat_scale;
            L = (lon[i] - rpcPara->long_off) / rpcPara->long_scale;
            H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;

            _ati_rpc_plh_coef(P, L, H, A);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -X);
            memset(a, 0, 59 * sizeof(double));
            memcpy(a, A, 20 * sizeof(double));
            memcpy(a + 40, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->b, A, 20);
            _ati_pnrml(a, 59, X, ata, atl, p*p);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -Y);
            memset(a, 0, 59 * sizeof(double));
            memcpy(a + 20, A, 20 * sizeof(double));
            memcpy(a + 40, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->d, A, 20);
            _ati_pnrml(a, 59, Y, ata, atl, p*p);
        }

        h = 1.0e-06;
        for (k = 0, n = 0; k<59; k++){
            ata[n] += h;
            n += (59 - k);
        }

        _ati_solve(ata, atl, x, 59, 59);

        maxv = fabs(rpcPara->a[0] - x[0]);
        h = fabs(rpcPara->c[0] - x[20]);
        if (h > maxv) maxv = h;
        for (i = 1; i<20; i++){
            h = fabs(rpcPara->a[i] - x[i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->c[i] - x[20 + i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->b[i] - x[40 + i - 1]);
            if (h > maxv) maxv = h;
        }

        rpcPara->b[0] = 1.0;
        rpcPara->d[0] = 1.0;
        memcpy(rpcPara->a, x, 20 * sizeof(double));
        memcpy(rpcPara->c, x + 20, 20 * sizeof(double));
        memcpy(rpcPara->b + 1, x + 40, 19 * sizeof(double));
        memcpy(rpcPara->d + 1, x + 40, 19 * sizeof(double));

        if (maxv < 1.0e-06 && iter > 0) break;
        iter++;
    } while (iter < 50);

    double dx, dy, ax, ay;
    sx = sy = ax = ay = 0.0;
    for (i = 0; i < num; i++){
        P = (lat[i] - rpcPara->lat_off) / rpcPara->lat_scale;
        L = (lon[i] - rpcPara->long_off) / rpcPara->long_scale;
        H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;
        _ati_rpc_plh_coef(P, L, H, A);

        X = _ati_rpc_scalar_product(rpcPara->a, A, 20) / _ati_rpc_scalar_product(rpcPara->b, A, 20);
        Y = _ati_rpc_scalar_product(rpcPara->c, A, 20) / _ati_rpc_scalar_product(rpcPara->d, A, 20);
        X = X*rpcPara->samp_scale + rpcPara->samp_off;
        Y = Y*rpcPara->line_scale + rpcPara->line_off;

        dx = samp[i] - X;
        dy = line[i] - Y;
        sx += (dx*dx);
        sy += (dy*dy);
        ax += dx;
        ay += dy;
    }
    sx = sqrt(sx / num);
    sy = sqrt(sy / num);
    ax /= num;
    ay /= num;

    printf("\nRPC Cal result\n");
    printf(" sx = %8.4lf   sy = %8.4lf\n", sx, sy);
    printf(" ax = %8.4lf   ay = %8.4lf\n", ax, ay);
};

// directly solve the 78 RPC parameters through the input 2D points and 3D points
static void _rpc_calc_rpc_78(const double* samp, const double* line, const  double* lat, const  double* lon, const  double* alt, int num, RpcPara* rpcPara, double& sx, double& sy, double &mx, double &my)
{
    memset(rpcPara, 0, sizeof(RpcPara));
    _ati_rpc_off_scale(samp, num, rpcPara->samp_off, rpcPara->samp_scale);
    _ati_rpc_off_scale(line, num, rpcPara->line_off, rpcPara->line_scale);
    _ati_rpc_off_scale(lat, num, rpcPara->lat_off, rpcPara->lat_scale);
    _ati_rpc_off_scale(lon, num, rpcPara->long_off, rpcPara->long_scale);
    _ati_rpc_off_scale(alt, num, rpcPara->height_off, rpcPara->height_scale);
    rpcPara->samp_scale /= 2.0;
    rpcPara->line_scale /= 2.0;
    rpcPara->lat_scale /= 2.0;
    rpcPara->long_scale /= 2.0;

    double a[78], ata[78 * 78], atl[78], x[78];
    double P, L, H, X, Y, A[20], B[20], p, h, maxv, minv = 999;
    RpcPara minvRpc; memset(&minvRpc, 0, sizeof(minvRpc));

    int i, iter = 0;
    do{
        memset(ata, 0, 78 * 78 * sizeof(double));
        memset(atl, 0, 78 * sizeof(double));
        memset(a, 0, 78 * sizeof(double));
        memset(x, 0, 78 * sizeof(double));

        for (i = 0; i < num; i++){
            X = (samp[i] - rpcPara->samp_off) / rpcPara->samp_scale;
            Y = (line[i] - rpcPara->line_off) / rpcPara->line_scale;
            P = (lat[i] - rpcPara->lat_off) / rpcPara->lat_scale;
            L = (lon[i] - rpcPara->long_off) / rpcPara->long_scale;
            H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;

            _ati_rpc_plh_coef(P, L, H, A);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -X);
            memset(a, 0, 78 * sizeof(double));
            memcpy(a, A, 20 * sizeof(double));
            memcpy(a + 20, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->b, A, 20);
            _ati_pnrml(a, 78, X, ata, atl, p*p);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -Y);
            memset(a, 0, 78 * sizeof(double));
            memcpy(a + 39, A, 20 * sizeof(double));
            memcpy(a + 59, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->d, A, 20);
            _ati_pnrml(a, 78, Y, ata, atl, p*p);
        }

        h = 1.0e-06;
        for (int k = 0, n = 0; k<78; k++){
            ata[n] += h;
            n += (78 - k);
        }

        _ati_solve(ata, atl, x, 78, 78);

        maxv = fabs(rpcPara->a[0] - x[0]);
        h = fabs(rpcPara->c[0] - x[39]);
        if (h > maxv) maxv = h;
        for (i = 1; i<20; i++){
            h = fabs(rpcPara->a[i] - x[i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->c[i] - x[39 + i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->b[i] - x[20 + i - 1]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->d[i] - x[59 + i - 1]);
            if (h > maxv) maxv = h;
        }

        rpcPara->b[0] = 1.0;
        rpcPara->d[0] = 1.0;
        memcpy(rpcPara->a, x, 20 * sizeof(double));
        memcpy(rpcPara->b + 1, x + 20, 19 * sizeof(double));
        memcpy(rpcPara->c, x + 39, 20 * sizeof(double));
        memcpy(rpcPara->d + 1, x + 59, 19 * sizeof(double));

        if (minv > maxv){ minv = maxv; memcpy(&minvRpc, rpcPara, sizeof(minvRpc)); }
        if (maxv < 1.0e-06 && iter > 0) break;
        iter++;
    } while (iter < 25);
    if (iter >= 25) memcpy(rpcPara, &minvRpc, sizeof(minvRpc));

    double dx, dy, ax, ay, mxx = 0, mxy = 0;
    sx = sy = ax = ay = 0.0;
    for (i = 0; i < num; i++){
        P = (lat[i] - rpcPara->lat_off) / rpcPara->lat_scale;
        L = (lon[i] - rpcPara->long_off) / rpcPara->long_scale;
        H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;
        _ati_rpc_plh_coef(P, L, H, A);

        X = _ati_rpc_scalar_product(rpcPara->a, A, 20) / _ati_rpc_scalar_product(rpcPara->b, A, 20);
        Y = _ati_rpc_scalar_product(rpcPara->c, A, 20) / _ati_rpc_scalar_product(rpcPara->d, A, 20);
        X = X*rpcPara->samp_scale + rpcPara->samp_off;
        Y = Y*rpcPara->line_scale + rpcPara->line_off;

        dx = samp[i] - X; if (fabs(dx) > fabs(mxx)) mxx = dx;
        dy = line[i] - Y; if (fabs(dy) > fabs(mxy)) mxy = dy;
        sx += (dx*dx);
        sy += (dy*dy);
        ax += dx;
        ay += dy;
    }
    sx = sqrt(sx / num); mx = mxx;
    sy = sqrt(sy / num); my = mxy;
    ax /= num;
    ay /= num;

    printf("\nRPC Cal result (iter=%d):\n", iter);
    printf("sgma: %.6lf \t  %.6lf\n", sx, sy);
    printf("avrg: %.6lf \t  %.6lf\n", ax, ay);
    printf("max : %.6lf \t  %.6lf\n", mxx, mxy);
};

// calculate 59 inverse RPC parameters through the input 2D points and 3D points
static void _rpc_calc_inv_rpc_59(double* samp, double* line, double* lat, double* lon, double* alt, int num, RpcPara* rpcPara, double& sx, double& sy)
{
    memset(rpcPara, 0, sizeof(RpcPara));
    _ati_rpc_off_scale(samp, num, rpcPara->samp_off, rpcPara->samp_scale);
    _ati_rpc_off_scale(line, num, rpcPara->line_off, rpcPara->line_scale);
    _ati_rpc_off_scale(lat, num, rpcPara->lat_off, rpcPara->lat_scale);
    _ati_rpc_off_scale(lon, num, rpcPara->long_off, rpcPara->long_scale);
    _ati_rpc_off_scale(alt, num, rpcPara->height_off, rpcPara->height_scale);
    rpcPara->samp_scale *= 2.0;
    rpcPara->line_scale *= 2.0;
    rpcPara->lat_scale *= 2.0;
    rpcPara->long_scale *= 2.0;

    double a[59], ata[59 * 59], atl[59], x[59];
    double P, L, H, X, Y, A[20], B[20], p, h, maxv;

    int i, n, k, iter = 0;
    do{
        memset(ata, 0, 59 * 59 * sizeof(double));
        memset(atl, 0, 59 * sizeof(double));
        memset(a, 0, 59 * sizeof(double));
        memset(x, 0, 59 * sizeof(double));

        for (i = 0; i < num; i++){
            X = (samp[i] - rpcPara->samp_off) / rpcPara->samp_scale;
            Y = (line[i] - rpcPara->line_off) / rpcPara->line_scale;
            P = (lat[i] - rpcPara->lat_off) / rpcPara->lat_scale;
            L = (lon[i] - rpcPara->long_off) / rpcPara->long_scale;
            H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;

            _ati_rpc_plh_coef(X, Y, H, A);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -L);
            memset(a, 0, 59 * sizeof(double));
            memcpy(a, A, 20 * sizeof(double));
            memcpy(a + 40, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->b, A, 20);
            _ati_pnrml(a, 59, L, ata, atl, p*p);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -P);
            memset(a, 0, 59 * sizeof(double));
            memcpy(a + 20, A, 20 * sizeof(double));
            memcpy(a + 40, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->d, A, 20);
            _ati_pnrml(a, 59, P, ata, atl, p*p);
        }

        h = 1.0e-06;
        for (k = 0, n = 0; k<59; k++){
            ata[n] += h;
            n += (59 - k);
        }

        _ati_solve(ata, atl, x, 59, 59);

        maxv = fabs(rpcPara->a[0] - x[0]);
        h = fabs(rpcPara->c[0] - x[20]);
        if (h > maxv) maxv = h;
        for (i = 1; i<20; i++){
            h = fabs(rpcPara->a[i] - x[i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->c[i] - x[20 + i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->b[i] - x[40 + i - 1]);
            if (h > maxv) maxv = h;
        }

        rpcPara->b[0] = 1.0;
        rpcPara->d[0] = 1.0;
        memcpy(rpcPara->a, x, 20 * sizeof(double));
        memcpy(rpcPara->c, x + 20, 20 * sizeof(double));
        memcpy(rpcPara->b + 1, x + 40, 19 * sizeof(double));
        memcpy(rpcPara->d + 1, x + 40, 19 * sizeof(double));

        if (maxv < 1.0e-06 && iter > 0) break;
        iter++;
    } while (iter < 50);

    double dx, dy, ax, ay;
    sx = sy = ax = ay = 0.0;
    for (i = 0; i < num; i++){
        X = (samp[i] - rpcPara->samp_off) / rpcPara->samp_scale;
        Y = (line[i] - rpcPara->line_off) / rpcPara->line_scale;
        H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;
        _ati_rpc_plh_coef(X, Y, H, A);

        L = _ati_rpc_scalar_product(rpcPara->a, A, 20) / _ati_rpc_scalar_product(rpcPara->b, A, 20);
        P = _ati_rpc_scalar_product(rpcPara->c, A, 20) / _ati_rpc_scalar_product(rpcPara->d, A, 20);
        L = L*rpcPara->long_scale + rpcPara->long_off;
        P = P*rpcPara->lat_scale + rpcPara->lat_off;

        dx = lon[i] - L;
        dy = lat[i] - P;
        sx += (dx*dx);
        sy += (dy*dy);
        ax += dx;
        ay += dy;
    }
    printf("\n");
    sx = sqrt(sx / num);
    sy = sqrt(sy / num);
    ax /= num;
    ay /= num;

    double l0 = rpcPara->lat_off*3.1415926535897932384626433832795 / 180.0;
    double ux = 6356752.314*3.1415926535897932384626433832795 / 180.0*cos(l0);
    double uy = 6378137.000*3.1415926535897932384626433832795 / 180.0;
    sx *= ux;
    sy *= uy;
    ax *= ux;
    ay *= uy;
    printf(" sx = %8.2lf   sy = %8.2lf\n", sx, sy);
    printf(" ax = %8.2lf   ay = %8.2lf\n", ax, ay);
};

// solve the 78 inverse RPC parameters through the input 2D points and 3D points
static void _rpc_calc_inv_rpc_78(const double* samp, const  double* line, const double* lat, const double* lon, const double* alt, int num, RpcPara* rpcPara, double& sx, double& sy)
{
    memset(rpcPara, 0, sizeof(RpcPara));
    _ati_rpc_off_scale(samp, num, rpcPara->samp_off, rpcPara->samp_scale);
    _ati_rpc_off_scale(line, num, rpcPara->line_off, rpcPara->line_scale);
    _ati_rpc_off_scale(lat, num, rpcPara->lat_off, rpcPara->lat_scale);
    _ati_rpc_off_scale(lon, num, rpcPara->long_off, rpcPara->long_scale);
    _ati_rpc_off_scale(alt, num, rpcPara->height_off, rpcPara->height_scale);
    rpcPara->samp_scale *= 2.0;
    rpcPara->line_scale *= 2.0;
    rpcPara->lat_scale *= 2.0;
    rpcPara->long_scale *= 2.0;

    double a[78], ata[78 * 78], atl[78], x[78];
    double P, L, H, X, Y, A[20], B[20], p, h, maxv;

    int i, k, n, iter = 0;
    do{
        memset(ata, 0, 78 * 78 * sizeof(double));
        memset(atl, 0, 78 * sizeof(double));
        memset(a, 0, 78 * sizeof(double));
        memset(x, 0, 78 * sizeof(double));

        for (i = 0; i < num; i++){
            X = (samp[i] - rpcPara->samp_off) / rpcPara->samp_scale;
            Y = (line[i] - rpcPara->line_off) / rpcPara->line_scale;
            P = (lat[i] - rpcPara->lat_off) / rpcPara->lat_scale;
            L = (lon[i] - rpcPara->long_off) / rpcPara->long_scale;
            H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;

            _ati_rpc_plh_coef(X, Y, H, A);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -L);
            memset(a, 0, 78 * sizeof(double));
            memcpy(a, A, 20 * sizeof(double));
            memcpy(a + 20, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->b, A, 20);
            _ati_pnrml(a, 78, L, ata, atl, p*p);

            memcpy(B, A, 20 * sizeof(double));
            _ati_rpc_numeric_product(B, 20, -P);
            memset(a, 0, 78 * sizeof(double));
            memcpy(a + 39, A, 20 * sizeof(double));
            memcpy(a + 59, B + 1, 19 * sizeof(double));
            p = (iter == 0) ? 1.0 : _ati_rpc_scalar_product(rpcPara->d, A, 20);
            _ati_pnrml(a, 78, P, ata, atl, p*p);
        }

        h = 1.0e-06;
        for (k = 0, n = 0; k<78; k++){
            ata[n] += h;
            n += (78 - k);
        }

        _ati_solve(ata, atl, x, 78, 78);

        maxv = fabs(rpcPara->a[0] - x[0]);
        h = fabs(rpcPara->c[0] - x[39]);
        if (h > maxv) maxv = h;
        for (i = 1; i<20; i++){
            h = fabs(rpcPara->a[i] - x[i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->c[i] - x[39 + i]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->b[i] - x[20 + i - 1]);
            if (h > maxv) maxv = h;
            h = fabs(rpcPara->d[i] - x[59 + i - 1]);
            if (h > maxv) maxv = h;
        }

        rpcPara->b[0] = 1.0;
        rpcPara->d[0] = 1.0;
        memcpy(rpcPara->a, x, 20 * sizeof(double));
        memcpy(rpcPara->b + 1, x + 20, 19 * sizeof(double));
        memcpy(rpcPara->c, x + 39, 20 * sizeof(double));
        memcpy(rpcPara->d + 1, x + 59, 19 * sizeof(double));

        if (maxv < 1.0e-06 && iter > 0) break;
        iter++;
    } while (iter < 50);

    double dx, dy, ax, ay;
    sx = sy = ax = ay = 0.0;
    for (i = 0; i < num; i++){
        X = (samp[i] - rpcPara->samp_off) / rpcPara->samp_scale;
        Y = (line[i] - rpcPara->line_off) / rpcPara->line_scale;
        H = (alt[i] - rpcPara->height_off) / rpcPara->height_scale;
        _ati_rpc_plh_coef(X, Y, H, A);

        L = _ati_rpc_scalar_product(rpcPara->a, A, 20) / _ati_rpc_scalar_product(rpcPara->b, A, 20);
        P = _ati_rpc_scalar_product(rpcPara->c, A, 20) / _ati_rpc_scalar_product(rpcPara->d, A, 20);
        L = L*rpcPara->long_scale + rpcPara->long_off;
        P = P*rpcPara->lat_scale + rpcPara->lat_off;

        dx = lon[i] - L;
        dy = lat[i] - P;
        sx += (dx*dx);
        sy += (dy*dy);
        ax += dx;
        ay += dy;
    }
    printf("\n");
    sx = sqrt(sx / num);
    sy = sqrt(sy / num);
    ax /= num;
    ay /= num;
    double l0 = rpcPara->lat_off*3.1415926535897932384626433832795 / 180.0;
    double ux = 6356752.314*3.1415926535897932384626433832795 / 180.0*cos(l0);
    double uy = 6378137.000*3.1415926535897932384626433832795 / 180.0;
    sx *= ux;
    sy *= uy;
    ax *= ux;
    ay *= uy;
    printf(" %8.2lf   %8.2lf\n", sx, sy);
    printf(" %8.2lf   %8.2lf\n", ax, ay);
};

////////////////////////////////////////////////////////////////////////////////////////////////
//CRpc

class CRpc
{
public:
    CRpc(){
        m_psTransform = NULL;
        memset(m_strDemPath, 0, sizeof(m_strDemPath));
    }
    ~CRpc(){
        if (m_psTransform) GDALDestroyRPCTransformer(m_psTransform);
    }
    void Reset()
    {
        ResetRpcInfo();
        if (m_psTransform) GDALDestroyRPCTransformer(m_psTransform); 
        m_psTransform = NULL;
    }
    bool AttachDemFile(const char* lpstrPathName)
    {
        if (!IsExist(lpstrPathName)) return false;
        strcpy(m_strDemPath, lpstrPathName);
        if (m_psTransform) GDALDestroyRPCTransformer(m_psTransform);
        m_psTransform = NULL;
        m_psTransform = CreateRpcTransformer(&m_rpcInfo, m_strDemPath);
        return true;
    }
    bool LoadPixFile(const char* lpstrPathName){
        GDALDataset* pImgSet = (GDALDataset *)GDALOpen(lpstrPathName, GA_ReadOnly);
        if (!pImgSet) {
            return false;
        }
        char **papszOptions = pImgSet->GetMetadata("RPC");
        int ret = GDALExtractRPCInfo(papszOptions, &m_rpcInfo);
        if (ret) {
            m_psTransform = CreateRpcTransformer(&m_rpcInfo, m_strDemPath);
        }
        GDALClose(pImgSet); //CSLDestroy(papszOptions);
        
        return true;
    }
    bool LoadRpcFile(const char* lpstrPathName)
    {
        Reset();
        if (strlen(lpstrPathName) < 3) return false;
        FILE *fRpc = fopen(lpstrPathName, "rt");    if (!fRpc) return false;

        int i; char chID[64]; double *rpc = (double*)(&m_rpcInfo);
        for (i = 0; i < 10; i++){ if (fscanf(fRpc, "%s%lf%s", chID, rpc + i, chID) != 3){ fclose(fRpc); return false; } }
        for (i = 0; i < 80; i++){ if (fscanf(fRpc, "%s%lf", chID, rpc + i + 10) != 2){ fclose(fRpc); return false; } }      
        fclose(fRpc);

        m_psTransform = CreateRpcTransformer(&m_rpcInfo,m_strDemPath);
        return true;
    }
    bool LoadRpbFile(const char* lpstrPathName)
    {
        Reset();
        if (strlen(lpstrPathName) < 3) return false;
        FILE *fRpc = fopen(lpstrPathName, "rt");    if (!fRpc) return false;
        int i; double *rpc = (double*)(&m_rpcInfo);
        // Read RPB ,PVL
        char strLine[256], strID[64], strEq[30], strVal[128];
        char *pGet = fgets(strLine, 256, fRpc);
        while (pGet){
            sscanf(strLine, "%s%s%s", strID, strEq, strVal);
            if (stricmp(strID, "BEGIN_GROUP") == 0 &&
                stricmp(strEq, "=") == 0 &&
                strnicmp(strVal, "IMAGE", 5) == 0){
                while (pGet){
                    pGet = fgets(strLine, 256, fRpc);
                    sscanf(strLine, "%s%s%s", strID, strEq, strVal);

                    if (stricmp(strID, "lineOffset") == 0) rpc[0] = atof(strVal); else
                    if (stricmp(strID, "sampOffset") == 0) rpc[1] = atof(strVal); else
                    if (stricmp(strID, "latOffset") == 0) rpc[2] = atof(strVal); else
                    if (stricmp(strID, "longOffset") == 0) rpc[3] = atof(strVal); else
                    if (stricmp(strID, "heightOffset") == 0) rpc[4] = atof(strVal); else
                    if (stricmp(strID, "lineScale") == 0) rpc[5] = atof(strVal); else
                    if (stricmp(strID, "sampScale") == 0) rpc[6] = atof(strVal); else
                    if (stricmp(strID, "latScale") == 0) rpc[7] = atof(strVal); else
                    if (stricmp(strID, "longScale") == 0) rpc[8] = atof(strVal); else
                    if (stricmp(strID, "heightScale") == 0) rpc[9] = atof(strVal); else
                    if (stricmp(strID, "lineNumCoef") == 0){
                        for (i = 10; (i < 30 && !feof(fRpc));){
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;
                            rpc[i++] = atof(strLine);
                        }
                    }
                    else
                    if (stricmp(strID, "lineDenCoef") == 0){
                        for (i = 30; (i < 50 && !feof(fRpc));){
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;
                            rpc[i++] = atof(strLine);
                        }
                    }
                    else
                    if (stricmp(strID, "sampNumCoef") == 0){
                        for (i = 50; (i < 70 && !feof(fRpc));){
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;
                            rpc[i++] = atof(strLine);
                        }
                    }
                    else
                    if (stricmp(strID, "sampDenCoef") == 0){
                        for (int i = 70; (i < 90 && !feof(fRpc));){
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;
                            rpc[i++] = atof(strLine);
                        }
                    }
                }
                break;
            }
            else
            if (stricmp(strID, "BEGIN_GROUP") == 0 &&
                stricmp(strEq, "=") == 0 &&
                strnicmp(strVal, "rationalFunctions", 17) == 0)
            {
                while (pGet)
                {
                    pGet = fgets(strLine, 256, fRpc);
                    sscanf(strLine, "%s%s%s", strID, strEq, strVal);

                    if (stricmp(strID, "lineOffset") == 0) rpc[0] = atof(strVal); else
                    if (stricmp(strID, "pixelOffset") == 0) rpc[1] = atof(strVal); else
                    if (stricmp(strID, "latitudeOffset") == 0) rpc[2] = atof(strVal); else
                    if (stricmp(strID, "longitudeOffset") == 0) rpc[3] = atof(strVal); else
                    if (stricmp(strID, "heightOffset") == 0) rpc[4] = atof(strVal); else
                    if (stricmp(strID, "lineScale") == 0) rpc[5] = atof(strVal); else
                    if (stricmp(strID, "pixelScale") == 0) rpc[6] = atof(strVal); else
                    if (stricmp(strID, "latitudeScale") == 0) rpc[7] = atof(strVal); else
                    if (stricmp(strID, "longitudeScale") == 0) rpc[8] = atof(strVal); else
                    if (stricmp(strID, "heightScale") == 0) rpc[9] = atof(strVal); else
                    if (stricmp(strID, "lineNumeratorCoefficients") == 0)
                    {
                        for (i = 10; (i < 30 && !feof(fRpc));)
                        {
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;

                            rpc[i++] = atof(strLine);
                        }
                    }
                    else
                    if (stricmp(strID, "lineDenominatorCoefficients") == 0)
                    {
                        for (i = 30; (i < 50 && !feof(fRpc));)
                        {
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;

                            rpc[i++] = atof(strLine);
                        }
                    }
                    else
                    if (stricmp(strID, "pixelNumeratorCoefficients") == 0)
                    {
                        for (i = 50; (i < 70 && !feof(fRpc));)
                        {
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;

                            rpc[i++] = atof(strLine);
                        }
                    }
                    else
                    if (stricmp(strID, "pixelDenominatorCoefficients") == 0)
                    {
                        for (int i = 70; (i < 90 && !feof(fRpc));)
                        {
                            pGet = fgets(strLine, 256, fRpc);
                            if (!pGet) continue; while (*pGet == ' ') pGet++;
                            if (strlen(pGet) < 3) continue;

                            rpc[i++] = atof(strLine);
                        }
                    }
                }
                break;
            }
            pGet = fgets(strLine, 256, fRpc);
        }

        fclose(fRpc);
        m_psTransform = CreateRpcTransformer(&m_rpcInfo, m_strDemPath);
        return true;
    }
    bool Load4File(const char* lpstrPathName)
    {
        if (strlen(lpstrPathName) < 3) return false;
        const char* pExt = strrchr(lpstrPathName, '.');
        if (pExt && (!strcmp(pExt + 1, "pix") || !strcmp(pExt + 1, "PIX"))){
            return LoadPixFile(lpstrPathName);
        }
        FILE *fRpc = fopen(lpstrPathName, "rt");    if (!fRpc) return false;

        int i; char chID[64];
        fscanf(fRpc, "%s", chID); fclose(fRpc);
        if (stricmp(chID, "LINE_OFF:") == 0){
            return LoadRpcFile(lpstrPathName);
        }
        else
        {
            return LoadRpbFile(lpstrPathName);
        }
        
        return false;
    }
    bool SaveRpcFile(const char* lpstrPathName)
    {
        if (strlen(lpstrPathName) < 3) return false;
        FILE *fAop = fopen(lpstrPathName, "wt");
        if (!fAop) return false;

        int i; double *rpc = (double*)&m_rpcInfo;
        fprintf(fAop, "LINE_OFF: %.3lf pixels\n", +rpc[0]);
        fprintf(fAop, "SAMP_OFF: %.3lf pixels\n", +rpc[1]);
        fprintf(fAop, "LAT_OFF: %.9lf degrees\n", +rpc[2]);
        fprintf(fAop, "LONG_OFF: %.9lf degrees\n", +rpc[3]);
        fprintf(fAop, "HEIGHT_OFF: %.6lf meters\n", +rpc[4]);
        fprintf(fAop, "LINE_SCALE: %.3lf pixels\n", +rpc[5]);
        fprintf(fAop, "SAMP_SCALE: %.3lf pixels\n", +rpc[6]);
        fprintf(fAop, "LAT_SCALE: %.9lf degrees\n", +rpc[7]);
        fprintf(fAop, "LONG_SCALE: %.9lf degrees\n", +rpc[8]);
        fprintf(fAop, "HEIGHT_SCALE: %lf meters\n", +rpc[9]);
        for (i = 1; i <= 20; i++)
            fprintf(fAop, "LINE_NUM_COEFF_%d: %.9E \n", i, +rpc[9 + i]);
        for (i = 1; i <= 20; i++)
            fprintf(fAop, "LINE_DEN_COEFF_%d: %.9E \n", i, +rpc[29 + i]);
        for (i = 1; i <= 20; i++)
            fprintf(fAop, "SAMP_NUM_COEFF_%d: %.9E \n", i, +rpc[49 + i]);
        for (i = 1; i <= 20; i++)
            fprintf(fAop, "SAMP_DEN_COEFF_%d: %.9E \n", i, +rpc[69 + i]);
        
        fclose(fAop);
        return true;
    }
    bool SaveRpbFile(const char* lpstrPathName)
    {
        if (strlen(lpstrPathName) < 3) return false;
        FILE *fp = fopen(lpstrPathName, "wt");
        if (!fp) return false;

        double *rpc = (double*)&m_rpcInfo;
        fprintf(fp, "satId = \"XXX\";\n");
        fprintf(fp, "bandId = \"XXX\";\n");
        fprintf(fp, "SpecId = \"XXX\";\n");
        fprintf(fp, "BEGIN_GROUP = IMAGE\n");
        fprintf(fp, "   errBias =   1.0;\n");
        fprintf(fp, "   errRand =   0.0;\n");
        fprintf(fp, "   lineOffset = %.15f;\n",     +rpc[0]);
        fprintf(fp, "   sampOffset = %.15f;\n",     +rpc[1]);
        fprintf(fp, "   latOffset = %.15f;\n",      +rpc[2]);
        fprintf(fp, "   longOffset = %.15f;\n",     +rpc[3]);
        fprintf(fp, "   heightOffset = %.15f;\n",   +rpc[4]);
        fprintf(fp, "   lineScale = %.15f;\n",      +rpc[5]);
        fprintf(fp, "   sampScale = %.15f;\n",      +rpc[6]);
        fprintf(fp, "   latScale = %.15f;\n",       +rpc[7]);
        fprintf(fp, "   longScale = %.15f;\n",      +rpc[8]);
        fprintf(fp, "   heightScale = %.15f;\n",    +rpc[9]);

        int i;
        fprintf(fp, "   lineNumCoef = (\n");
        for ( i = 0; i < 19; ++i)
            fprintf(fp, "\t\t\t%.15e,\n", +rpc[10 + i]);
        fprintf(fp, "\t\t\t%.15e);\n", + rpc[10 + i]);

        fprintf(fp, "   lineDenCoef = (\n");
        for ( i = 0; i < 19; ++i)
            fprintf(fp, "\t\t\t%.15e,\n", +rpc[30 + i]);
        fprintf(fp, "\t\t\t%.15e);\n", +rpc[30 + i]);

        fprintf(fp, "   sampNumCoef = (\n");
        for ( i = 0; i < 19; ++i)
            fprintf(fp, "\t\t\t%.15e,\n", +rpc[50 + i]);
        fprintf(fp, "\t\t\t%.15e);\n", +rpc[50 + i]);

        fprintf(fp, "   sampDenCoef = (\n");
        for ( i = 0; i < 19; ++i)
            fprintf(fp, "\t\t\t%.15e,\n", +rpc[70 + i]);
        fprintf(fp, "\t\t\t%.15e);\n", +rpc[70 + i]);

        fprintf(fp, "END_GROUP = IMAGE\nEND;\n");

        fclose(fp);
        return true;
    }
    void Calc_rpc_78(const double* samp, const  double* line, const  double* lat, const double* lon, const  double* alt, int num, double *sx0, double *sy0, double *mx0, double *my0)
    {
        double sx = 0, sy = 0, mx = 0, my = 0;
        _rpc_calc_rpc_78(samp, line, lat, lon, alt, num, (RpcPara*)&m_rpcInfo, sx, sy, mx, my);
        if (sx0) *sx0 = sx; if (sy0) *sy0 = sy;
        if (mx0) *mx0 = mx; if (my0) *my0 = my;
    }
    void GetRpcPara(RpcPara* par)
    {
        memcpy(par,&m_rpcInfo,sizeof(RpcPara));
    }
    //if dem is loaded, alt is relative height to dem 
    void RPCGrd2Pho_GDAL(double* lat, double* lon, double* alt,int nPtNum)
    {
        int* panSuccess = new int[nPtNum];
        GDALRPCTransform(m_psTransform,1,nPtNum,lon,lat,alt,panSuccess);
        delete[] panSuccess;
    }
    void RPCPhoZ2Grd_GDAL(double* px, double* py, double* alt, int nPtNum)
    {
        int* panSuccess = new int[nPtNum];
        GDALRPCTransform(m_psTransform, 0, nPtNum, px, py, alt, panSuccess);
        delete[] panSuccess;
    }
    void RPCGrd2Pho_WHU(double lat,  double lon,  double alt, double *px, double *py)
    {
        _rpc_geoedtic_to_pxy((RpcPara*)&m_rpcInfo, lat, lon, alt, px, py);
    }

    void RPCPhoZ2Grd_WHU(double px,  double py,  double gz, double *lat, double *lon, double *alt)
    {
        _rpc_pxy_Z_to_geoedtic((RpcPara*)&m_rpcInfo, px, py, gz, lat, lon, alt);
    }

    const char* GetLastErrorMsg() { return CPLGetLastErrorMsg(); }
protected:
    void ResetRpcInfo()
    {
        memset(&m_rpcInfo, 0, sizeof(double)*90);
        m_rpcInfo.dfMIN_LONG = -180.0;
        m_rpcInfo.dfMIN_LAT = -90.0;
        m_rpcInfo.dfMAX_LONG = 180.0;
        m_rpcInfo.dfMAX_LAT = 90.0;
    }
    
    static void* CreateRpcTransformer(GDALRPCInfo* psRpcInfo, const char* lpstrDemPath)
    {
        char **papszOptions = NULL;
        if(strlen(lpstrDemPath)>3)  papszOptions = CSLSetNameValue(papszOptions, "RPC_DEM", lpstrDemPath);

        return GDALCreateRPCTransformer(psRpcInfo, 0, 0.005, papszOptions);
    }
private:
    GDALRPCInfo m_rpcInfo;
    void*       m_psTransform;
    char        m_strDemPath[512];
};

CRpcBase::CRpcBase()
{
    m_rpc = new CRpc;
    m_pDxyGrid = NULL;
    memset(m_aop6, 0, sizeof(m_aop6));
}


CRpcBase::~CRpcBase()
{
    if (m_rpc) delete m_rpc; m_rpc = NULL;
    if (m_pDxyGrid) delete m_pDxyGrid; m_pDxyGrid = NULL;
}

void CRpcBase::Reset()
{
    memset(m_aop6, 0, sizeof(m_aop6));
    if (m_pDxyGrid) delete m_pDxyGrid; m_pDxyGrid = NULL;
    m_gridC = m_gridR = 0;
    m_gridDx = m_gridDy = 0;
}

bool CRpcBase::Load4File(const char* lpstrPathName)
{
    Reset();
    if (m_rpc->Load4File(lpstrPathName)){
        FILE* fRpc = fopen(lpstrPathName, "rt");    if (!fRpc) return false;

        char strLine[256];  char chID[64];
        char *pGet = fgets(strLine, 256, fRpc);
        while (pGet)
        {
            sscanf(strLine, "%s", chID);
            if (strnicmp(chID, "RFM_CORRECTION_PARAMETERS", 25) == 0){ for (int i = 0; i < 6; i++) fscanf(fRpc, "%s%lf\n", chID, m_aop6 + i); } 
            if (feof(fRpc)) break; 
            pGet = fgets(strLine, 256, fRpc);
        }
        fclose(fRpc);

        char strGrd[512]; strcpy(strGrd, lpstrPathName); strcat(strGrd, ".grd");
        FILE *fGrd = fopen(strGrd, "rb");
        if (fGrd){
            fread(&m_gridC, sizeof(int), 1, fGrd);
            fread(&m_gridR, sizeof(int), 1, fGrd);
            fread(&m_gridDx, sizeof(int), 1, fGrd);
            fread(&m_gridDy, sizeof(int), 1, fGrd);
            if (NewGrid(m_gridC, m_gridR, m_gridDx, m_gridDy)){
                fread(m_pDxyGrid, m_gridC*m_gridR*sizeof(int), 1, fGrd);
            }
            else{ m_gridC = m_gridR = 0; }
            fclose(fGrd);
        }
        return true;
    }
    return false;
}

bool CRpcBase::Save2File(const char* lpstrPathName, double *sx /* = NULL */, double *sy /* = NULL */, double *mx /* = NULL */, double *my /* = NULL */)
{
    if (strlen(lpstrPathName)<3) return false;
    char strExt[10];    strcpy(strExt, strrchr(lpstrPathName, '.') + 1); strlwr(strExt);
    bool ret = false;
    if (!strcmp(strExt, "rpb")) ret = m_rpc->SaveRpbFile(lpstrPathName);
    else ret = m_rpc->SaveRpcFile(lpstrPathName);

    if (!ret){
        return false;
    }
    
    FILE* fAop = fopen(lpstrPathName, "a+");
    if (fAop){
        double *ab6 = m_aop6;
        fprintf(fAop, "RFM_CORRECTION_PARAMETERS:\n");
        fprintf(fAop, "R0:  %.9E\n", ab6[0]);
        fprintf(fAop, "R1:  %.9E\n", ab6[1]);
        fprintf(fAop, "R2:  %.9E\n", ab6[2]);
        fprintf(fAop, "C0:  %.9E\n", ab6[3]);
        fprintf(fAop, "C1:  %.9E\n", ab6[4]);
        fprintf(fAop, "C2:  %.9E\n", ab6[5]);
        fclose(fAop);
    }
    if (m_gridC > 0 && m_gridR > 0){
        char strGrd[512]; strcpy(strGrd, lpstrPathName); strcat(strGrd, ".grd");
        FILE *fGrd = fopen(strGrd, "wb");
        if (fGrd){
            fwrite(&m_gridC, sizeof(int), 1, fGrd);
            fwrite(&m_gridR, sizeof(int), 1, fGrd);
            fwrite(&m_gridDx, sizeof(int), 1, fGrd);
            fwrite(&m_gridDy, sizeof(int), 1, fGrd);
            fwrite(m_pDxyGrid, m_gridC*m_gridR*sizeof(int), 1, fGrd);
            fclose(fGrd);
        }
    }

    if (sx&&sy&&mx&&my){
        char strRes[512]; strcpy(strRes, lpstrPathName); strcat(strRes, ".res");
        FILE *fRes = fopen(strRes, "wt");
        if (fRes){
            fprintf(fRes, "RPC Cal Result:\n");
            fprintf(fRes, "sgma: %.6lf \t  %.6lf\n", *sx, *sy);
            fprintf(fRes, "max : %.6lf \t  %.6lf\n", *mx, *my);
            fclose(fRes);
        }
    }
    return true;
}

void CRpcBase::Calc_rpc_78(const double* samp, const  double* line, const  double* lat, const double* lon, const  double* alt, int num, double *sx /* = NULL */, double *sy /* = NULL */, double *mx /* = NULL */, double *my /* = NULL */)
{
    memset(m_aop6, 0, sizeof(m_aop6));
    m_rpc->Calc_rpc_78(samp, line, lat, lon, alt, num, sx, sy, mx, my);
}

void CRpcBase::RPCGrd2Pho(double lat, double lon, double alt, double *px, double *py)
{
    m_rpc->RPCGrd2Pho_WHU(lat, lon, alt, px, py);
    if (m_pDxyGrid){
        double dx = 0, dy = 0; GetDxy(*px, *py, &dx, &dy);
        *px += dx; *py += dy;
    }
}
void CRpcBase::RPCGrd2Pho(const double* lat, const double* lon, double* alt, int nPtNum, double *px, double *py)
{
    memcpy(px, lon, sizeof(double)*nPtNum);
    memcpy(py, lat, sizeof(double)*nPtNum);
    m_rpc->RPCGrd2Pho_GDAL(py,px,alt,nPtNum);
    if (m_pDxyGrid){
        double dx = 0, dy = 0; 
        for (int i = 0; i < nPtNum;i++){
            dx = 0, dy = 0;
            GetDxy(*px, *py, &dx, &dy);
            *px += dx; *py += dy;
            px++;   py++;
        }
        
    }
}
void CRpcBase::RPCPhoZ2Grd(double px, double py, double gz, double *lat, double *lon, double *alt)
{
    if (m_pDxyGrid){
        double dx = 0, dy = 0; GetDxy(px, py, &dx, &dy);
        px -= dx; py -= dx;
    }
    m_rpc->RPCPhoZ2Grd_WHU(px, py, gz, lat, lon, alt);
}
void CRpcBase::RPCPhoZ2Grd(const double* px, const  double* py, const  double* gz, int nPtNum, double *lat, double *lon, double *alt)
{
    RPCPhoZ2Grd(px, py, gz, nPtNum, lat, lon, alt, false);
}

void CRpcBase::RPCPhoZ2Grd(const double* x, const double* y, const double* gz, int nPtNum, double *lat, double *lon, double *alt, bool bAop6)
{
    memcpy(lon, x, sizeof(double)*nPtNum);
    memcpy(lat, y, sizeof(double)*nPtNum);
    memcpy(alt, gz, sizeof(double)*nPtNum);
    
    if (bAop6){
        double *pLon = lon, *plat = lat;
        for (int i = 0; i < nPtNum; i++) RPC_ixy_to_pxy(*pLon, *plat, pLon, plat);
    }
    if (m_pDxyGrid){
        double *pLon = lon, *plat = lat;
        double dx = 0, dy = 0;
        for (int i = 0; i < nPtNum; i++){
            dx = 0, dy = 0;
            GetDxy(*pLon, *plat, &dx, &dy);
            *pLon -= dx; *plat -= dy;
            pLon++; plat++;
        }
    }

    m_rpc->RPCPhoZ2Grd_GDAL(lon, lat, alt, nPtNum);
}

void CRpcBase::GetRpcPara(RpcPara* para)
{
    m_rpc->GetRpcPara(para);
}

bool CRpcBase::AttachDemFile(const char* lpstrPathName)
{
    return m_rpc->AttachDemFile(lpstrPathName);
}

void RPCPho2Grd(CRpcBase *rpcL, double pxl, double pyl,
    CRpcBase *rpcR, double pxr, double pyr,
    double *lat, double *lon, double *alt){
    RpcPara paraL, paraR;
    rpcL->GetRpcPara(&paraL);   rpcR->GetRpcPara(&paraR);
    _rpc_pxy_to_geoedtic(&paraL, pxl, pyl, &paraR, pxr, pyr, lat, lon, alt);
}
