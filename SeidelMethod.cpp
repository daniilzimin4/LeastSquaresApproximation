// Daniil Zimin
// d.zimin@innopolis.university
// tg: daniilzimin4

#include <bits/stdc++.h>

#include <utility>

using namespace std;

typedef long long ll;
typedef long double ld;

mt19937 rnd(time(0));

#define ioss ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
#define file freopen("input.txt", "r", stdin); freopen("output.txt", "w", stdout);
#define pb push_back
#define all(x) x.begin(),x.end()
#define rall(x) x.rbegin(),x.rend()

template<class T>
ll upmax(T &a, T b) { return (b > a) ? a = b, 1 : 0; }

template<class T>
ll upmin(T &a, T b) { return (b < a) ? a = b, 1 : 0; }

const ll inf = 1e18 + 7;
const ld EPS = 1e-6;

class Matrix {
public:
    int row{}, column{};
    vector<vector<ld>> matrix;

    Matrix() {
        row = 0, column = 0;
    }

    Matrix(int n, int m) {
        row = n;
        column = m;
        matrix.resize(n, vector<ld>(m, 0));
    }

    friend istream &operator>>(istream &input, Matrix &cur_matrix) {
        for (int i = 0; i < cur_matrix.row; i++)
            for (int j = 0; j < cur_matrix.column; j++)
                input >> cur_matrix.matrix[i][j];
        return input;
    }

    friend ostream &operator<<(ostream &output, Matrix &cur_matrix) {
        for (int i = 0; i < cur_matrix.row; i++) {
            for (int j = 0; j < cur_matrix.column; j++)
                output << fixed << setprecision(4) << (abs(cur_matrix.matrix[i][j]) < EPS ? 0 : cur_matrix.matrix[i][j])
                       << " ";
            output << endl;
        }
        return output;
    }

    Matrix &operator=(const Matrix &cur_matrix) {
        column = cur_matrix.column;
        row = cur_matrix.row;
        for (int i = 0; i < column; i++)
            for (int j = 0; j < row; j++)
                matrix[i][j] = cur_matrix.matrix[i][j];
        return *this;
    }

    Matrix operator+(const Matrix &second) {
        if (this->row != second.row || this->column != second.column) {
            cout << "Error: the dimensional problem occurred\n";
            exit(0);
        }
        Matrix ans(this->row, this->column);
        for (int i = 0; i < this->row; i++)
            for (int j = 0; j < second.column; j++)
                ans.matrix[i][j] = this->matrix[i][j] + second.matrix[i][j];
        return ans;
    }

    Matrix operator-(const Matrix &second) {
        if (this->row != second.row || this->column != second.column) {
            cout << "Error: the dimensional problem occurred\n";
            exit(0);
        }
        Matrix ans(this->row, this->column);
        for (int i = 0; i < this->row; i++)
            for (int j = 0; j < second.column; j++)
                ans.matrix[i][j] = this->matrix[i][j] - second.matrix[i][j];
        return ans;
    }

    Matrix operator*(const Matrix &second) {
        if (this->column != second.row) {
            cout << "Error: the dimensional problem occurred\n";
            exit(0);
        }
        Matrix ans(this->row, second.column);
        for (int i = 0; i < this->row; i++)
            for (int j = 0; j < second.column; j++)
                for (int l = 0; l < this->matrix[i].size(); l++)
                    ans.matrix[i][j] += this->matrix[i][l] * second.matrix[l][j];
        return ans;
    }

    Matrix transpose() {
        Matrix ans(this->column, this->row);
        for (int i = 0; i < this->row; i++)
            for (int j = 0; j < this->column; j++)
                ans.matrix[j][i] = this->matrix[i][j];
        return ans;
    }
};

class SquareMatrix : public Matrix {
public:
    explicit SquareMatrix(int x) {
        row = column = x;
        matrix.resize(x, vector<ld>(x));
    }

    explicit SquareMatrix(const Matrix &ss) {
        if (ss.row != ss.column) {
            cout << "Error: the dimensional problem occurred\n";
            exit(0);
        }
        row = column = ss.row;
        matrix.resize(row, vector<ld>(row));
        for (int i = 0; i < row; i += 1) {
            for (int j = 0; j < column; j += 1) {
                matrix[i][j] = ss.matrix[i][j];
            }
        }
    }

    SquareMatrix() {
        column = row = 0;
    }

    friend istream &operator>>(istream &input, SquareMatrix &o) {
        for (int i = 0; i < o.row; i += 1)
            for (int j = 0; j < o.row; j += 1)
                input >> o.matrix[i][j];
        return input;
    }

    friend ostream &operator<<(ostream &output, SquareMatrix &o) {
        for (int i = 0; i < o.row; i += 1) {
            for (int j = 0; j < o.row; j += 1)
                output << fixed << setprecision(4) << (abs(o.matrix[i][j]) < EPS ? 0 : o.matrix[i][j]) << ' ';
            output << endl;
        }
        return output;
    }

    SquareMatrix operator+(SquareMatrix second) {
        Matrix ans = (*((Matrix *) this)) + (*((Matrix *) &second));
        return *(SquareMatrix *) &ans;
    }

    SquareMatrix operator-(SquareMatrix second) {
        Matrix ans = (*((Matrix *) this)) - (*((Matrix *) &second));
        return *(SquareMatrix *) &ans;
    }

    SquareMatrix operator*(SquareMatrix second) {
        Matrix ans = (*((Matrix *) this)) * (*((Matrix *) &second));
        return *(SquareMatrix *) &ans;
    }

    SquareMatrix transpose() {
        Matrix res = ((Matrix) *this).transpose();
        return *(SquareMatrix *) &res;
    }


    int getMaxInColumn(int col) {
        int ans = col;
        for (int j = col; j < row; j += 1)
            ans = (abs(matrix[ans][col]) < abs(matrix[j][col])) ? j : ans;
        return ans;
    }
};

class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix() {
        row = column = 0;
    }

    explicit IdentityMatrix(int n) {
        row = column = n;
        matrix.resize(n, vector<ld>(n, 0.0));
        for (int i = 0; i < n; i++)
            matrix[i][i] = 1.0;
    }
};

class PermutationMatrix : public IdentityMatrix {
public:
    PermutationMatrix(const SquareMatrix &x, int i, int j) : IdentityMatrix(x.column) {
        swap(matrix[i - 1], matrix[j - 1]);
    }
};

class ColumnVector : public Matrix {
public:
    ColumnVector() = default;

    explicit ColumnVector(int n) {
        row = n;
        column = 1;
        matrix.resize(row, vector<ld>(column));
    }

    friend istream &operator>>(istream &input, ColumnVector &o) {
        for (int i = 0; i < o.row; i += 1)
            for (int j = 0; j < o.column; j += 1)
                input >> o.matrix[i][j];
        return input;
    }

    friend ostream &operator<<(ostream &output, ColumnVector &o) {
        for (int i = 0; i < o.row; i += 1) {
            for (int j = 0; j < o.column; j += 1)
                output << fixed << setprecision(4)
                       << (abs(o.matrix[i][j]) < EPS ? 0 : o.matrix[i][j]);
            output << endl;
        }
        return output;
    }
};

ld getDeterminant(SquareMatrix a) {
    int n = a.column;
    int k = 1;
    auto coutStep = [&k, &a](int type) {
        if (type == 1) {
            cout << "step #" << k++ << ": permutation\n";
        } else {
            cout << "step #" << k++ << ": elimination\n";
        }
        cout << a;
    };
    for (int l = 0; l < n - 1; l += 1) {
        // Find pivots
        int tmp = a.getMaxInColumn(l);
        if (tmp != l) {
            PermutationMatrix p(a, 1 + tmp, 1 + l);
            a = p * a;
            coutStep(1);
        }
        // Convert to REF
        for (int i = l + 1; i < n; i += 1) {
            if (abs(a.matrix[i][l]) > EPS) {
                ld delta = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j += 1) {
                    a.matrix[i][j] = a.matrix[i][j] - delta * a.matrix[l][j];
                }
                coutStep(2);
            }
        }
    }
    ld ans = 1;
    for (int i = 0; i < a.row; i += 1) {
        ans = ans * a.matrix[i][i];
    }
    cout << "result:\n";
    return ans;
}

void solveEquation(SquareMatrix a, ColumnVector b, int comments = 0) {
    int n = a.column, k = 0;
    auto printStep = [&a, &b, &k, n, comments](int type) {
        if (!comments)return;
        if (type == 0) {
            cout << "step #" << k++ << ":\n";
        } else if (type == 1) {
            cout << "step #" << k++ << ": elimination\n";
        } else if (type == 2) {
            cout << "step #" << k++ << ": permutation\n";
        } else if (type == 3) {
            cout << "Diagonal normalization:\n";
        }
        cout << a;
        cout << b;
    };
    printStep(0);
    for (int l = 0; l < n - 1; l += 1) {
        // Find pivots
        int tmp = a.getMaxInColumn(l);
        if (l != tmp) {
            PermutationMatrix p(a, 1 + tmp, 1 + l);
            a = p * a;
            swap(b.matrix[tmp], b.matrix[l]);
            printStep(2);
        }
        // Convert to REF
        for (int i = l + 1; i < n; i += 1) {
            if (abs(a.matrix[i][l]) > EPS) {
                ld delta = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j += 1) {
                    a.matrix[i][j] = a.matrix[i][j] - delta * a.matrix[l][j];
                }
                b.matrix[i][0] = b.matrix[i][0] - delta * b.matrix[l][0];
                printStep(1);
            }
        }
    }
    // Check for INF answer or NO (if REF consists of all zeroes in last row and last element in column)
    int ok = 0;
    for (int i = 0; i < n; i += 1)
        ok = ok || (abs(a.matrix[n - 1][i]) > EPS);
    if (!ok) {
        if (abs(b.matrix[n - 1][0]) < EPS) {
            cout << "INF";
        } else {
            cout << "NO";
        }
        exit(0);
    }
    // Convert to Identity matrix
    for (int l = n - 1; l >= 0; l -= 1) {
        for (int i = l - 1; i >= 0; i -= 1) {
            if (abs(a.matrix[i][l]) > EPS) {
                ld delta = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j += 1) {
                    a.matrix[i][j] = a.matrix[i][j] - delta * a.matrix[l][j];
                }
                b.matrix[i][0] = b.matrix[i][0] - delta * b.matrix[l][0];
                printStep(1);
            }
        }
    }
    // Normalization column vector
    for (int i = 0; i < n; i += 1) {
        b.matrix[i][0] = b.matrix[i][0] / a.matrix[i][i];
        a.matrix[i][i] = 1.0;
    }
    printStep(3);
    if (comments)cout << "result:\n";
    cout << b;
}

void InversePrint(SquareMatrix a) {
    int n = a.column, k = 0;
    SquareMatrix res = IdentityMatrix(n);
    auto printStep = [&a, &res, &k, n](int type) {
        if (type == 0) {
            cout << "step #" << k++ << ": Augmented Matrix\n";
        } else if (type == 1) {
            cout << "step #" << k++ << ": permutation\n";
        } else if (type == 2) {
            cout << "step #" << k++ << ": elimination\n";
        } else if (type == 3) {
            cout << "Diagonal normalization:\n";
        }
        for (int i = 0; i < n; i += 1) {
            cout << fixed << setprecision(2);
            for (int j = 0; j < n; j += 1) {
                cout << (abs(a.matrix[i][j]) < EPS ? 0 : a.matrix[i][j]) << ' ';
            }
            for (int j = 0; j < n; j += 1) {
                cout << (abs(res.matrix[i][j]) < EPS ? 0 : res.matrix[i][j]) << ' ';
            }
            cout << '\n';
        }
    };
    printStep(0);
    cout << "Direct way:\n";
    for (int l = 0; l < n - 1; l += 1) {
        // Find pivots
        int tmp = a.getMaxInColumn(l);
        if (tmp != l) {
            PermutationMatrix p(a, 1 + tmp, 1 + l);
            a = p * a, res = p * res;
            printStep(1);
        }
        // Convert to REF
        for (int i = l + 1; i < n; i += 1) {
            if (abs(a.matrix[i][l]) > EPS) {
                ld delta = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j += 1) {
                    a.matrix[i][j] = a.matrix[i][j] - delta * a.matrix[l][j];
                    res.matrix[i][j] = res.matrix[i][j] - delta * res.matrix[l][j];
                }
                printStep(2);
            }
        }
    }
    // Convert to Identity matrix
    cout << "Way back:\n";
    for (int l = n - 1; l >= 0; l -= 1) {
        for (int i = l - 1; i >= 0; i -= 1) {
            if (abs(a.matrix[i][l]) > EPS) {
                ld delta = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j += 1) {
                    a.matrix[i][j] = a.matrix[i][j] - delta * a.matrix[l][j];
                    res.matrix[i][j] = res.matrix[i][j] - delta * res.matrix[l][j];
                }
                printStep(2);
            }
        }
    }
    // Normalization identity matrix
    for (int i = 0; i < n; i += 1) {
        for (int j = 0; j < n; j += 1) {
            res.matrix[i][j] = res.matrix[i][j] / a.matrix[i][i];
        }
        a.matrix[i][i] = 1.0;
    }
    printStep(3);
    cout << "result:\n";
    cout << res;
}

SquareMatrix Inverse(SquareMatrix a) {
    int n = a.column, k = 0;
    SquareMatrix res = IdentityMatrix(n);
    for (int l = 0; l < n - 1; l += 1) {
        // Find pivots
        int tmp = a.getMaxInColumn(l);
        if (tmp != l) {
            PermutationMatrix p(a, 1 + tmp, 1 + l);
            a = p * a, res = p * res;
        }
        // Convert to REF
        for (int i = l + 1; i < n; i += 1) {
            if (abs(a.matrix[i][l]) > EPS) {
                ld delta = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j += 1) {
                    a.matrix[i][j] = a.matrix[i][j] - delta * a.matrix[l][j];
                    res.matrix[i][j] = res.matrix[i][j] - delta * res.matrix[l][j];
                }
            }
        }
    }
    // Convert to Identity matrix
    for (int l = n - 1; l >= 0; l -= 1) {
        for (int i = l - 1; i >= 0; i -= 1) {
            if (abs(a.matrix[i][l]) > EPS) {
                ld delta = a.matrix[i][l] / a.matrix[l][l];
                for (int j = 0; j < n; j += 1) {
                    a.matrix[i][j] = a.matrix[i][j] - delta * a.matrix[l][j];
                    res.matrix[i][j] = res.matrix[i][j] - delta * res.matrix[l][j];
                }
            }
        }
    }
    // Normalization identity matrix
    for (int i = 0; i < n; i += 1) {
        for (int j = 0; j < n; j += 1) {
            res.matrix[i][j] = res.matrix[i][j] / a.matrix[i][i];
        }
        a.matrix[i][i] = 1.0;
    }
    return res;
}

pair<SquareMatrix, SquareMatrix> LU(const SquareMatrix& a) {
    int n = a.row;
    SquareMatrix B(n), C = a;
    for (int i = 0; i < n; i += 1) {
        for (int j = i; j < n; j += 1) {
            B.matrix[j][i] = C.matrix[j][i] / C.matrix[i][i];
        }
        for (int j = i + 1; j < n; j += 1) {
            for (int k = i; k < n; k += 1) {
                C.matrix[j][k] -= C.matrix[i][k] * B.matrix[j][i];
            }
        }
    }
    return {B, C};
}

void solve() {
    int n;
    cin >> n;
    SquareMatrix A(n);
    cin >> A;
    int m;
    cin >> m;
    ColumnVector b(m);
    cin >> b;
    ld eps;
    cin >> eps;
    SquareMatrix D(n);
    for (int i = 0; i < n; i += 1) {
        D.matrix[i][i] = A.matrix[i][i];
    }
    for (int i = 0; i < n; i += 1) {
        ld sum = 0;
        for (int j = 0; j < n; j += 1) {
            sum += A.matrix[i][j];
        }
        if (abs(A.matrix[i][i]) < abs(sum - A.matrix[i][i])) {
            cout << "The method is not applicable!\n";
            exit(0);
        }
    }


    SquareMatrix D_1 = Inverse(D);

    SquareMatrix alpha = IdentityMatrix(n) - (D_1 * A);
    SquareMatrix B1(n), C1(n);
    for (int i = 0; i < n; i += 1) {
        for (int j = 0; j < n; j += 1) {
            if (j < i) {
                B1.matrix[i][j] = alpha.matrix[i][j];
            }
            if (j > i) {
                C1.matrix[i][j] = alpha.matrix[i][j];
            }
        }
    }
    SquareMatrix B(n), C(n);
    for(int i = 0; i < n; i += 1){
        for(int j = 0; j < n; j += 1){
            if(j <= i)B.matrix[i][j] = A.matrix[i][j];
            else C.matrix[i][j] = A.matrix[i][j];
        }
    }
    SquareMatrix B_1 = Inverse(B);
    Matrix beta = ((Matrix) D_1 * b);
    cout << "beta:\n" << beta;
    cout << "alpha:\n" << alpha;
    cout << "B:\n" << B1;
    cout << "C:\n" << C1;
    SquareMatrix I_B = IdentityMatrix(n) - B1;
    cout << "I-B:\n" << I_B;
    SquareMatrix I_B_1 = Inverse(I_B);
    cout << "(I-B)_-1:\n" << I_B_1;

    Matrix X = beta;
    cout << "x(" << 0 << "):\n" << X;
    ld epstmp = inf;
    int id = 1;
    while (epstmp > eps) {
        Matrix Xnew = (Matrix) B_1 * (b - (Matrix) C * X);
        Matrix tmpX = Xnew - X;
        epstmp = 0;
        for (int j = 0; j < m; j += 1) {
            epstmp += (tmpX.matrix[j][0]) * (tmpX.matrix[j][0]);
            X.matrix[j][0] = Xnew.matrix[j][0];
        }
        epstmp = sqrt(epstmp);
        cout << "e: " << epstmp << '\n';
        cout << "x(" << id++ << "):\n" << Xnew;
    }
}

int main() {
//    ioss
//    file
//    (a += b) -> return value
    int t = 1;
//    cin >> t;
    while (t--) {
        solve();
    }
    // printf("Time taken: %.10fs\n", (double)(clock() - CLOCK_START)/CLOCKS_PER_SEC);
    return 0;
}
