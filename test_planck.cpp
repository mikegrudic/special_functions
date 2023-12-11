#include "planck.c"
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

int main()
{
    double x1 = 1.;
    double x2 = x1;
    std::ofstream myfile;
    myfile.open("x_vs_planck.dat");
    while (x2 < 100.)
    {
        myfile << std::setprecision(16) << std::scientific << x2 << " " << planck_integral(x1, x2) << "\n";
        x2 *= 1.1;
    }
    myfile << std::setprecision(16) << std::scientific << INFINITY << " " << planck_integral(x1, INFINITY) << "\n";
    myfile.close();
    return 0;
}
