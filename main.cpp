#include <cassert>
#include <complex>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

typedef complex<double> cd;
typedef Triplet<cd> Tcd;

const int L = 10;
const int n = 2;
const double J = 1.0;

const int N = pow(n, L);
const double s = (n - 1.0) / 2.0;

MatrixXcd generate_single(vector<Tcd> S)
{
    SparseMatrix<cd> M(n, n);
    M.setFromTriplets(S.begin(), S.end());
    return M;
}

SparseMatrix<cd> generate(int i, vector<Tcd> S)
{
    assert(i < L);
    SparseMatrix<cd> M(N, N);
    int i_shift = pow(n, i);
    int upper_shift = i_shift * n;
    int lower_size = i_shift;
    int upper_size = N / upper_shift;
    for (Tcd x : S)
    {
        assert(x.row() < n);
        assert(x.col() < n);
        int row = x.row() * i_shift;
        int col = x.col() * i_shift;
        for (int iu = 0; iu < upper_size; iu++)
        {
            int upper = iu * upper_shift;
            int upper_row = upper + row;
            int upper_col = upper + col;
            for (int il = 0; il < lower_size; il++)
            {
                M.insert(upper_row + il, upper_col + il) = x.value();
            }
        }
    }
    return M;
}

int main()
{
    clock_t start = clock();

    vector<Tcd> I, Sx, Sy, Sz;
    for (int i = 0; i < n; i++)
    {
        I.push_back(Tcd(i, i, 1.0));
        double sz = s - i;
        Sz.push_back(Tcd(i, i, sz));
        if (i >= 1)
        {
            double a = sqrt(i * (n - i)) / 2.0;
            Sx.push_back(Tcd(i - 1, i, a));
            Sx.push_back(Tcd(i, i - 1, a));
            Sy.push_back(Tcd(i - 1, i, a * cd(0.0, -1.0)));
            Sy.push_back(Tcd(i, i - 1, a * cd(0.0, 1.0)));
        }
    }

    SparseMatrix<cd> H(N, N);
    for (int i = 0; i < L; i++)
    {
        int j = i + 1;
        if (j >= L)
        {
            j %= L;
        }
        SparseMatrix<cd> SiSj = generate(i, Sx) * generate(j, Sx) + generate(i, Sy) * generate(j, Sy) + generate(i, Sz) * generate(j, Sz);
        H += SiSj;
    }
    H *= J;

    SparseMatrix<cd> Sx_total(N, N), Sy_total(N, N), Sz_total(N, N);
    for (int i = 0; i < L; i++)
    {
        Sx_total += generate(i, Sx);
        Sy_total += generate(i, Sy);
        Sz_total += generate(i, Sz);
    }
    SparseMatrix<cd> S_total = Sx_total * Sx_total + Sy_total * Sy_total + Sz_total * Sz_total;

    cout << fixed;

    // SelfAdjointEigenSolver<SparseMatrix<cd>> es_s;
    // es_s.compute(10 * S_total + Sz_total);
    // cout << "eigenvalues(S, Sz):" << endl;
    // cout << es_s.eigenvalues()(seqN(0, 20)) << endl;
    // cout << "eigenvectors:" << endl;
    // cout << es_s.eigenvectors() << endl;

    SelfAdjointEigenSolver<SparseMatrix<cd>> es;
    es.compute(H);
    clock_t end = clock();
    std::cout << "calculation time: " << (double)(end - start) / CLOCKS_PER_SEC << "s." << endl;

    cout << "eigenvalues:" << endl;
    cout << es.eigenvalues()(seqN(0, 8)) << endl;
    cout << "spin gap: " << es.eigenvalues()[1] - es.eigenvalues()[0] << endl;
    // cout << "eigenvectors:" << endl;
    // cout << (es_s.eigenvectors().adjoint() * es.eigenvectors())(all, seqN(0, 4)) << endl;

    return 0;
}
