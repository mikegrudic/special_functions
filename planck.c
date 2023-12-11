#include "math.h"

double planck_lower(double x);
double planck_upper(double x);
double planck_integral(double x1, double x2);

double planck_upper(double x)
{
    // Series solution for upper Planck integral from x to infinity evaluated to machine
    // precision. Most efficient for large x.
    double tol = 2.220446049250313e-16;
    double series_sum = 0;
    double err = 1e100;
    double expx_inv = exp(-x);
    double expx_inv_power = expx_inv; // will cumulatively multiply at each iteration
    int n = 1;
    while (fabs(err) > tol * series_sum)
    {
        double nx = n * x;
        err = (6 + nx * (6 + nx * (3 + nx))) * expx_inv_power / (n * n * n * n);
        expx_inv_power *= expx_inv;
        series_sum += err;
        n += 1;
    }
    return series_sum * 0.15398973382026503;
}

double planck_lower(double x)
{
    // Series solution for Planck function integral from 0 to x,
    // most efficient for large x, valid to machine precision on [0,3]
    double coeffs[] = {
        0.3333333333333333,
        -0.125,
        0.016666666666666666,
        -0.0001984126984126984,
        3.6743092298647855e-6,
        -7.515632515632516e-8,
        1.6059043836821615e-9,
        -3.522793425791662e-11,
        7.872080312167458e-13,
        -1.784042261222412e-14,
        4.088600979179926e-16,
        -9.455950863295921e-18,
        2.203601131344092e-19,
        -5.168320254004638e-21,
        1.2188644964239545e-22,
        -2.888231428076628e-24,
        6.87258318890207e-26,
        -1.641368762534915e-27,
        3.932898582742878e-29,
        -9.451269078629001e-31,
        2.2772522578280595e-32,
        -5.500052129536349e-34,
        1.3312603916626966e-35,
        -3.22862741376232e-37,
        7.844404337661609e-39,
        -1.909088837773861e-40,
    };
    int num_coeffs = 25;

    double x_sqr = x * x;
    double val = coeffs[num_coeffs - 1] * x_sqr;

    int n = num_coeffs - 1;
    while (n >= 3)
    {
        val = x_sqr * (val + coeffs[n]);
        n--;
    }
    while (n >= 0)
    {
        val = x * (val + coeffs[n]);
        n--;
    }
    val *= x_sqr;

    return val * 0.15398973382026503;
}

double planck_integral(double x1, double x2)
{
    // Returns the definite integral of the normalized (frequency) Planck
    // function f = norm * x^3/(exp(x)-1) from x1 to x2, accurate to machine
    // precision

    double sign;
    if (x2 < x1)
    {
        // assume x1 < x2 throughout, and just remember the original parity
        double temp = x1;
        x1 = x2;
        x2 = temp;
        sign = -1;
    }
    else
    {
        sign = 1;
    }
    const double cutoff_upper_lower = 3; // where to switch between small x and large x series approx

    // initialize flag for if we need to take care with precision
    int upper_precision = 0;

    double f1, f2;
    if (x1 == 0)
    {
        f1 = 0;
    }
    else
    {
        if (x1 < cutoff_upper_lower)
        {
            f1 = planck_lower(x1);
        }
        else
        {
            upper_precision = 1;
            f1 = planck_upper(x1);
        }
    }

    if (x2 == INFINITY)
    {
        if (upper_precision)
        {
            return f1 * sign;
        }
        return (1 - f1) * sign;
    }
    if (x2 < cutoff_upper_lower)
    {
        f2 = planck_lower(x2);
        return (f2 - f1) * sign;
    }
    f2 = planck_upper(x2);
    if (upper_precision)
    {
        return (f1 - f2) * sign;
    }
    return (1 - f1 - f2) * sign;
}