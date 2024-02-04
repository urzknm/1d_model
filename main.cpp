#include <cassert>
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

typedef Triplet<double> Td;

int digit_count(int m, int n);
vector<int> index_2d_search_lower(int dl, int n, vector<int> &nup_count_dl);
vector<int> index_2d_search_upper(int du, int nup, int n, vector<int> &nup_count_dl);
int index_2d_search(int i, int n, vector<int> &index_l, vector<int> &index_u);

vector<Td> Lz(int n);
vector<Td> Lp(int n);
vector<Td> Lm(int n);

const double J = 1.0;

int L;

MatrixXcd generate_single(vector<Td> S, int n)
{
    SparseMatrix<double> M(n, n);
    M.setFromTriplets(S.begin(), S.end());
    return M;
}

SparseMatrix<double> generate(int i, vector<Td> S, int n)
{
    assert(i < L);
    int N = pow(n, L);
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

// Σ_<i,j> S1_i * S2_j, (S1, S2) = (Sz, Sz), (Sp, Sm), (Sm, Sp)
SparseMatrix<double> generate_couple_sum(vector<Td> S1, vector<Td> S2, int dc_all, int n, vector<int> &index_l, vector<int> &index_u)
{
    int N_all = pow(n, L);
    int N_rest = pow(n, L - 2);
    int index_max = index_u.back() + index_l.back() + 1;
    SparseMatrix<double> M(index_max, index_max);
    for (Td x : S1)
    {
        for (Td y : S2)
        {
            assert(x.row() < n);
            assert(x.col() < n);
            assert(y.row() < n);
            assert(y.col() < n);
            int row = (n * x.row() + y.row());
            int col = (n * x.col() + y.col());
            double val = (x.value() * y.value());
            int dc_ij = x.row() + y.row(); // digit_count(row, n);
            assert(dc_ij == x.col() + y.col());
            row *= N_rest;
            col *= N_rest;
            for (int ir = 0; ir < N_rest; ir++)
            {
                if (dc_ij + digit_count(ir, n) == dc_all)
                {
                    int row_all = row + ir;
                    int col_all = col + ir;
                    for (int shift = 0; shift < L; shift++)
                    {
                        M.coeffRef(index_2d_search(row_all, n, index_l, index_u), index_2d_search(col_all, n, index_l, index_u)) += val;
                        row_all = ((row_all % n) * N_all + row_all) / n;
                        col_all = ((col_all % n) * N_all + col_all) / n;
                    }
                }
            }
        }
    }
    return M;
}

int main()
{
    clock_t start = clock();

    int n;
    cout << "input (n L): ";
    cin >> n >> L;

    assert(pow(n, L) < INT_MAX);
    int N = pow(n, L);
    double s = (n - 1.0) / 2.0;

    vector<Td> Sz = Lz(n);
    vector<Td> Sp = Lp(n);
    vector<Td> Sm = Lm(n);

    int sz_all = 0;
    int dc_all = sz_all + s * L;

    vector<int> dc_l;
    vector<int> index_l = index_2d_search_lower(L / 2, n, dc_l);
    vector<int> index_u = index_2d_search_upper(L - (L / 2), dc_all, n, dc_l);

    SparseMatrix<double> H = generate_couple_sum(Sz, Sz, dc_all, n, index_l, index_u) + generate_couple_sum(Sp, Sm, dc_all, n, index_l, index_u) + generate_couple_sum(Sm, Sp, dc_all, n, index_l, index_u);
    // cout << H;
    /*
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
    */
    H *= J;

    /*
    SparseMatrix<double> Sp_total(N, N), Sm_total(N, N), Sz_total(N, N);
    for (int i = 0; i < L; i++)
    {
        Sp_total += generate(i, Sp, N);
        Sm_total += generate(i, Sm, N);
        Sz_total += generate(i, Sz, N);
    }
    SparseMatrix<double> S_total = Sp_total * Sm_total + Sm_total * Sp_total + Sz_total * Sz_total;
    */

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

    cout << "eigenvalues on Sz=" << sz_all << ":" << endl;
    cout << es.eigenvalues()(seqN(0, min((int)(H.rows()), 8))) << endl;
    assert(n > 1);
    assert(L >= 1);
    cout << "spin gap: Δ = " << es.eigenvalues()[1] - es.eigenvalues()[0] << "." << endl;
    // cout << "eigenvectors:" << endl;
    // cout << (es_s.eigenvectors().adjoint() * es.eigenvectors())(all, seqN(0, 4)) << endl;

    return 0;
}
