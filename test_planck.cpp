#include "planck.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

int main()
{
    double x1 = 0.1;
    double x2 = x1;
    int p = 3;

    std::ofstream myfile;
    myfile.open("x_vs_planck.dat");
    while (x2 < 100.)
    {
        myfile << std::setprecision(16) << std::scientific << x1 << " " << x2 << " " << planck_integral(x1, x2, p) << "\n";
        x2 *= 1.01;
    }
    myfile << std::setprecision(16) << std::scientific << x1 << " " << INFINITY << " " << planck_integral(x1, INFINITY, p) << "\n";
    myfile.close();
    return 0;
}
