
#include <cstdio>
#include <assert.h>
#include <stdlib.h>
#include <sys/time.h>

// A simple Matrix of unsigned integers laid out row-major in a 1D array.
// M is number of rows and N is number of columns.
struct Matrix
{
    int* data = nullptr;
    int M = 0;
    int N = 0;
};

// A simple function using for displaying a matrix
void ShowMatrix(const Matrix& m) {
    for (int row = 0; row < m.M; row++) {
        for (int col = 0; col < m.N; col++) {
            printf("%6d ", m.data[row * m.N + col]);
        }
        printf("\n");
    }
}

// Create and return a new matrix MxN, initialized with running integers
// starting with gap and advancing by gap to the next element. For example,
// CreateNewMatrix(2, 3, 5) will produce the matrix:
//
//   5  10  15
//  20  15  30
//
Matrix newMatrix(int M, int N, int gap) {
    Matrix matrix;
    matrix.M = M;
    matrix.N = N;
    matrix.data = new int[M * N];

    int value = gap;
    for (int i = 0; i < M * N; i++) {
        matrix.data[i] = value;
        value += gap;
    }
    return matrix;
}

/*ADDITION*/

// Add matrix x into matrix y (y += x) and store the result under the 3rd parameter
__attribute__((noinline))
__attribute__((optimize("no-tree-vectorize")))
void AddMatrixByRow(const Matrix& y, const Matrix& x, const Matrix &result) {
    assert(y.M == x.M);
    assert(y.N == x.N);

    for (int row = 0; row < y.M; row++) {
        for (int col = 0; col < y.N; col++) {
            result.data[row * y.N + col] += x.data[row * x.N + col];
        }
    }
}

// Add matrix x into matrix y (y += x) and store the result under the 3rd parameter
__attribute__((noinline))
__attribute__((optimize("no-tree-vectorize")))
void AddMatrixByCol(const Matrix& y, const Matrix& x, const Matrix &result) {
    assert(y.M == x.M);
    assert(y.N == x.N);

    for (int col = 0; col < y.N; col++) {
        for (int row = 0; row < y.M; row++) {
            result.data[row * y.N + col] += x.data[row * x.N + col];
        }
    }
}

/*SUBTRACTION*/

// Sub matrix x into matrix y (y -= x) and store the result under the 3rd parameter
__attribute__((noinline))
__attribute__((optimize("no-tree-vectorize")))
void SubMatrixByRow(const Matrix& y, const Matrix& x, const Matrix &result) {
    assert(y.M == x.M);
    assert(y.N == x.N);

    for (int row = 0; row < y.M; row++) {
        for (int col = 0; col < y.N; col++) {
            result.data[row * y.N + col] += x.data[row * x.N + col];
        }
    }
}

// Sub matrix x into matrix y (y -= x) and store the result under the 3rd parameter
__attribute__((noinline))
__attribute__((optimize("no-tree-vectorize")))
void SubMatrixByCol(const Matrix& y, const Matrix& x, const Matrix &result) {
    assert(y.M == x.M);
    assert(y.N == x.N);

    for (int col = 0; col < y.N; col++) {
        for (int row = 0; row < y.M; row++) {
            result.data[row * y.N + col] -= x.data[row * x.N + col];
        }
    }
}

/*MULTIPLICATION*/

// Mult matrix x into matrix y (y *= x) and store the result under the 3rd parameter
__attribute__((noinline))
__attribute__((optimize("no-tree-vectorize")))
void MulMatrixByRow(const Matrix& y, const Matrix& x, const Matrix &result) {
    assert(y.M == x.M);
    assert(y.N == x.N);

    for(int row = 0; row < y.M; row++) {
        for(int col = 0; col < x.N; col++) {
            result.data[row * x.N + col] = 0;
            for(int k = 0; k < x.M; k++) {
                result.data[row * x.N + col] = result.data[row * x.N + col] + y.data[row * x.N + k] * x.data[k * x.N + col];
            }
        }
    }
}

// Mult matrix x into matrix y (y *= x) and store the result under the 3rd parameter
__attribute__((noinline))
__attribute__((optimize("no-tree-vectorize")))
void MulMatrixByCol(const Matrix& x, const Matrix& y, const Matrix &result) {
    assert(y.M == x.M);
    assert(y.N == x.N);

    for (int col = 0; col < x.N; col++) {
        for (int row = 0; row < y.M; row++) {
            for(int k = 0; k < x.M; k++) {
                result.data[row * x.N + col] = result.data[row * x.N + col] + y.data[row * x.N + k] * x.data[k * x.N + col];
            }
        }
    }
}


// Return the current time
static unsigned long long get_timestamp ()
{
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (unsigned long long)now.tv_sec * 1000000;
}

// return the time needed for making : (m1 * m2) && (m1 + m2) && (m1 - m2)
// by using rowMajor calculation
double benchRowMajor(const Matrix& m1, const Matrix& m2, const Matrix& m3) {
    unsigned long long t_start = get_timestamp();
    AddMatrixByRow(m1, m2, m3);
    SubMatrixByRow(m1, m2, m3);
    MulMatrixByRow(m1, m2, m3);
    unsigned long long t_end = get_timestamp();

    return ((t_start - t_end) / 1000000.0L);
}

// return the time needed for making : (m1 * m2) && (m1 + m2) && (m1 - m2)
// by using columnMajor calculation
double benchColumnsMajor(const Matrix& m1, const Matrix& m2, const Matrix& m3) {
    unsigned long long t_start = get_timestamp();
    AddMatrixByCol(m1, m2, m3);
    SubMatrixByCol(m1, m2, m3);
    MulMatrixByCol(m1, m2, m3);
    unsigned long long t_end = get_timestamp();

    return ((t_start - t_end) / 1000000.0L);
}

/*MAIN*/
int main(int argc, char **agv) {

    Matrix m1 = newMatrix(500, 500, 5);
    Matrix m2 = newMatrix(500, 500, 3);
    Matrix m3 = newMatrix(500, 500, 0);

    double rowMajorExecutionTime = benchRowMajor(m1, m2, m3);
    double columnMajorExecutionTime = benchColumnsMajor(m1, m2, m3);
    double difference = rowMajorExecutionTime - columnMajorExecutionTime;

    printf("row-major execution time : %fs\n", rowMajorExecutionTime);
    printf("column-major execution time :%fs\n", columnMajorExecutionTime);
    printf("difference : %fs\n", difference < 0 ? -difference : difference);
}
