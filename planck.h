// Functions to evaluate definite integrals of the normalized, dimensionless
// frequency Planck function f = x^p/(exp(x)-1)
#include "math.h"
#include <iostream>

using namespace std;
double planck_integral(double x1, double x2, int p);
double planck_gaussquad(double x1, double x2, int p);
double planck_upper_series(double x1, double x2, int p);

// Returns the definite integral of the normalized (frequency) Planck
// function f = x^p/(exp(x)-1) from x1 to x2, accurate to machine
// precision for any valid floating point arguments on [0,inf)
double planck_integral(double x1, double x2, int p) {
    double sign = 1.;
    if (x1 > x2) { // assume x2 is the larger in magnitude
        double tmp = x2;
        x2 = x1;
        x1 = tmp;
        sign = -1;
    }
    const double cutoff = 3.0;
    if (x1 > cutoff) {
        return sign * planck_upper_series(x1, x2, p);
    }
    if (x2 < cutoff) {
        return sign * planck_gaussquad(x1, x2, p);
    }
    return sign * (planck_gaussquad(x1, cutoff, p) +
                   planck_upper_series(cutoff, x2, p));
}

// Series solution for Planck function integral from 0 to x,
// most efficient for large x, accurate to machine precision on [0,3],
// inaccurate at larger values
double planck3_lower(double x) {
    double coeffs[] = {
        0.3333333333333333,     -0.125,
        0.016666666666666666,   -0.0001984126984126984,
        3.6743092298647855e-6,  -7.515632515632516e-8,
        1.6059043836821615e-9,  -3.522793425791662e-11,
        7.872080312167458e-13,  -1.784042261222412e-14,
        4.088600979179926e-16,  -9.455950863295921e-18,
        2.203601131344092e-19,  -5.168320254004638e-21,
        1.2188644964239545e-22, -2.888231428076628e-24,
        6.87258318890207e-26,   -1.641368762534915e-27,
        3.932898582742878e-29,  -9.451269078629001e-31,
        2.2772522578280595e-32, -5.500052129536349e-34,
        1.3312603916626966e-35, -3.22862741376232e-37,
        7.844404337661609e-39,  -1.909088837773861e-40,
    };
    int num_coeffs = 25;

    double x_sqr = x * x;
    double val = coeffs[num_coeffs - 1] * x_sqr;
    int n = num_coeffs - 1;
    while (n >= 3) {
        val = x_sqr * (val + coeffs[n]);
        n--;
    }
    while (n >= 0) {
        val = x * (val + coeffs[n]);
        n--;
    }
    val *= x_sqr;
    return val;
}

double planck_gaussquad(double x1, double x2, int p) {
    if (x1 == x2) {
        return 0;
    }
    const int NUM_POINTS = 12;
    double roots[] = {
        -0.9815606342467192, -0.9041172563704749, -0.7699026741943047,
        -0.5873179542866175, -0.3678314989981801, -0.125233408511469,
        0.125233408511469,   0.3678314989981801,  0.5873179542866175,
        0.7699026741943047,  0.9041172563704749,  0.9815606342467192,
    };

    double weights[] = {
        0.0471753363865132, 0.1069393259953178, 0.1600783285433461,
        0.2031674267230657, 0.2334925365383547, 0.2491470458134026,
        0.2491470458134026, 0.2334925365383547, 0.2031674267230657,
        0.1600783285433461, 0.1069393259953178, 0.0471753363865132,
    };
    double integral = 0.0, mu = 2;

    for (int i = 0; i < NUM_POINTS; i++) {
        double x = 0.5 * (x2 - x1) * (roots[i] + 1.0) + x1;
        double xpow = 1.0;
        for (int k = 0; k < p; k++) {
            xpow *= x;
        }
        integral += xpow / expm1(x) * weights[i];
    }
    return integral * (x2 - x1) / mu;
}

double planck_upper_series(double x1, double x2, int p) {
    /* Series solution for Planck integral x^p/(exp(x)-1) from x1 to
    x2 evaluated to machine precision. Most efficient for large x.

    Parameters
    ----------
    x1: float
        Lower integral limit
    x2: float, optional
        Upper integration limit (default: inf)
    p: int, optional
        Exponent of Planck function (default: 3)

    Returns
    -------
    Planck integral evaluated to machine precision
    */
    const double tol = 2.220446049250313e-16;
    double series_sum = 0.0, term = 1e100, term1 = 1.0, term2 = 1.0, sign = 1.0;

    if (x1 > x2) { // assume x2 is the larger in magnitude
        double tmp = x2;
        x2 = x1;
        x1 = tmp;
        sign = -1;
    }
    // will cumulatively multiply exp(-x) at each iteration
    double expx1_inv = exp(-x1), expx1_inv_power = expx1_inv;
    double expx2_inv = 0.0, expx2_inv_power = 0.0;
    if (x2 < INFINITY) {
        expx2_inv = expx2_inv_power = exp(-x2);
    }

    // now sum the terms. Note this is the general case; we could simplify
    // if we know p a priori and make this a little faster
    double nx1, nx2, nprod;
    int m, n = 1;
    while (fabs(term) > tol * fabs(series_sum)) {
        if (x2 == INFINITY) {
            nx1 = n * x1;
            term1 = 1.0;
            m = 1;
            nprod = n;
            for (int i = p; i > 0; i--) {
                m *= i;
                term1 = term1 * nx1 + m;
                nprod *= n;
            }
            term1 /= nprod;
            term1 *= expx1_inv_power;
            expx1_inv_power *= expx1_inv;
            term2 = 0;
        } else { // both x1 and x2 are finite
            nx1 = n * x1;
            nx2 = n * x2;
            term1 = 1.0, term2 = 1.0;
            m = 1;
            nprod = n;
            for (int i = p; i > 0; i--) {
                m *= i;
                term1 = term1 * nx1 + m;
                term2 = term2 * nx2 + m;
                nprod *= n;
            }
            term1 /= nprod;
            term2 /= nprod;
            term1 *= expx1_inv_power;
            term2 *= expx2_inv_power;
            expx1_inv_power *= expx1_inv;
            expx2_inv_power *= expx2_inv;
        }
        term = term1 - term2;
        series_sum += term;
        n += 1;
    }
    return sign * series_sum;
}