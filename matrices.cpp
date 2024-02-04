#include <cassert>
#include <cmath>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

typedef Triplet<double> Td;

vector<Td> Lz(int n)
{
    double s = (n - 1.0) / 2.0;
    vector<Td> Lz;
    for (int i = 0; i < n; i++)
    {
        double sz = s - i;
        Lz.push_back(Td(i, i, sz));
    }
    return Lz;
}

vector<Td> Lp(int n)
{
    double s = (n - 1.0) / 2.0;
    vector<Td> Lp;
    for (int i = 1; i < n; i++)
    {
        double a = sqrt(i * (n - i) / 2.0);
        Lp.push_back(Td(i - 1, i, a));
    }
    return Lp;
}

vector<Td> Lm(int n)
{
    double s = (n - 1.0) / 2.0;
    vector<Td> Lm;
    for (int i = 1; i < n; i++)
    {
        double a = sqrt(i * (n - i) / 2.0);
        Lm.push_back(Td(i, i - 1, a));
    }
    return Lm;
}
