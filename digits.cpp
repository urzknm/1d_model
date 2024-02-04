#include <cassert>
#include <cmath>
#include <vector>

using namespace std;

// 自然数mのn進数表記での桁の数字の和
int digit_count(int m, int n)
{
    assert(m >= 0);
    int count = 0;
    while (m != 0)
    {
        count += m % n;
        m /= n;
    }
    return count;
}

/*
// digit_count[d][k] d桁以下のn進数で桁の数字の和がkであるものは何通りあるかのリスト
vector<vector<int>> digit_count(int d_max, int n)
{
    vector<vector<int>> count(0, vector<int>(1, 1));
    for (int d = 1; d <= d_max; d++)
    {
        vector<int> count_d((n - 1) * d + 1, 0);
        for (int top_num = 0; top_num < n; top_num++)
        {
            for (int k = top_num; k - top_num < count[d - 1].size(); k++)
            {
                count_d[k] += count[d - 1][k - top_num];
            }
        }
        count.push_back(count_d);
    }
    return count;
}
*/

vector<int> index_2d_search_lower(int dl, int n, vector<int> &digit_count_dl)
{
    // vector<int> digit_count_dl((n - 1) * dl + 1, 0);
    digit_count_dl.assign((n - 1) * dl + 1, 0);
    int Nl = pow(n, dl);
    vector<int> list(Nl);
    for (int il = 0; il < Nl; il++)
    {
        int dc_il = digit_count(il, n);
        list[il] = digit_count_dl[dc_il];
        digit_count_dl[dc_il]++;
    }
    return list;
}

vector<int> index_2d_search_upper(int du, int dc_all, int n, vector<int> &digit_count_dl)
{
    int Nu = pow(n, du);
    vector<int> list(Nu);
    int ju = 0;
    for (int iu = 0; iu < Nu; iu++)
    {
        list[iu] = ju;
        int dc_iu = digit_count(iu, n);
        int dc_il = dc_all - dc_iu;
        if (0 <= dc_il && dc_il < digit_count_dl.size())
            ju += digit_count_dl[dc_il];
    }
    return list;
}

int index_2d_search(int i, int n, vector<int> &index_l, vector<int> &index_u)
{
    int Nl = index_l.size();
    int il = i % Nl;
    int iu = i / Nl;
    assert(iu < index_u.size());
    return index_u[iu] + index_l[il];
}
