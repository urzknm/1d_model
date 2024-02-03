#include <cassert>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

typedef Triplet<double> Td;

const double J = 1.0;

int n;
int L;

MatrixXcd generate_single(vector<Td> S)
{
    SparseMatrix<double> M(n, n);
    M.setFromTriplets(S.begin(), S.end());
    return M;
}

SparseMatrix<double> generate(int i, vector<Td> S, int N)
{
    assert(i < L);
    // n進数で右からi桁目の値が、i番目の格子に対応
    SparseMatrix<double> M(N, N);
    int i_shift = pow(n, i);
    int upper_shift = i_shift * n;
    int lower_size = i_shift;
    int upper_size = N / upper_shift;
    for (Td x : S)
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

    cout << "input (n L): ";
    cin >> n >> L;

    assert(pow(n, L) < INT_MAX);
    int N = pow(n, L);
    double s = (n - 1.0) / 2.0;

    vector<Td> I, Sz, Sp, Sm;
    for (int i = 0; i < n; i++)
    {
        I.push_back(Td(i, i, 1.0));
        double sz = s - i;
        Sz.push_back(Td(i, i, sz));
        if (i >= 1)
        {
            double a = sqrt(i * (n - i) / 2.0);
            Sp.push_back(Td(i - 1, i, a));
            Sm.push_back(Td(i, i - 1, a));
        }
    }

    SparseMatrix<double> H(N, N);
    for (int i = 0; i < L; i++)
    {
        int j = i + 1;
        if (j >= L)
        {
            j %= L;
        }
        SparseMatrix<double> SiSj = generate(i, Sp, N) * generate(j, Sm, N) + generate(i, Sm, N) * generate(j, Sp, N) + generate(i, Sz, N) * generate(j, Sz, N);
        H += SiSj;
    }
    H *= J;

    SparseMatrix<double> Sp_total(N, N), Sm_total(N, N), Sz_total(N, N);
    for (int i = 0; i < L; i++)
    {
        Sp_total += generate(i, Sp, N);
        Sm_total += generate(i, Sm, N);
        Sz_total += generate(i, Sz, N);
    }
    SparseMatrix<double> S_total = Sp_total * Sm_total + Sm_total * Sp_total + Sz_total * Sz_total;

    cout << fixed;

    // SelfAdjointEigenSolver<SparseMatrix<double>> es_s;
    // es_s.compute(10 * S_total + Sz_total);
    // cout << "eigenvalues(S, Sz):" << endl;
    // cout << es_s.eigenvalues()(seqN(0, 20)) << endl;
    // cout << "eigenvectors:" << endl;
    // cout << es_s.eigenvectors() << endl;

    SelfAdjointEigenSolver<SparseMatrix<double>> es;
    es.compute(H);
    clock_t end = clock();
    std::cout << "calculation time: " << (double)(end - start) / CLOCKS_PER_SEC << "s." << endl;

    cout << "eigenvalues:" << endl;
    cout << es.eigenvalues()(seqN(0, min(N, 8))) << endl;
    cout << "spin gap: Δ = " << es.eigenvalues()[1] - es.eigenvalues()[0] << "." << endl;
    // cout << "eigenvectors:" << endl;
    // cout << (es_s.eigenvectors().adjoint() * es.eigenvectors())(all, seqN(0, 4)) << endl;

    return 0;
}
