#include "SIMDMatrixes.h"
#include <iostream>
#include <iomanip>
#include <array>

using namespace std;

void SIMDMatrixes::set_input_elements_of_massiv()
{
    while (true)
    {
        cout << "Введите m: " << endl;
        cin >> m;
        cout << "Введите p: " << endl;
        cin >> p;
        cout << "Введите q: " << endl;
        cin >> q;
        cout << "Введите количество проц-ых элементов: " << endl;
        cin >> process_element;
        if (m == 0 || p == 0 || q == 0 || process_element == 0) {
            cout << "Ошибка!Повторите ввод данных\n";
        }
        else break;

    }
}

double SIMDMatrixes::get_random_elements()
{
    double number = (double)(rand() % 20001) / 10000 - 1;
    return number;

}

void SIMDMatrixes::filling_matrixs()
{
    A = new double* [p];
    for (int i = 0; i < p; i++)
        A[i] = new double[m];

    for (int i = 0; i < p; i++)
        for (int j = 0; j < m; j++)
            A[i][j] = get_random_elements();

    B = new double* [m];
    for (int i = 0; i < m; i++)
        B[i] = new double[q];

    for (int i = 0; i < m; i++)
        for (int j = 0; j < q; j++)
            B[i][j] = get_random_elements();

    E = new double* [1];
    for (int i = 0; i < 1; i++)
        E[i] = new double[m];

    for (int i = 0; i < 1; i++)
        for (int j = 0; j < m; j++)
            E[i][j] = get_random_elements();

    G = new double* [p];
    for (int i = 0; i < p; i++)
        G[i] = new double[q];

    for (int i = 0; i < p; i++)
        for (int j = 0; j < q; j++)
            G[i][j] = get_random_elements();

    F = new double** [p];
    for (int i = 0; i < p; i++)
        F[i] = new double* [q];
    for (int i = 0; i < p; i++)
        for (int j = 0; j < q; j++)
            F[i][j] = new double[m];

    D = new double** [p];
    for (int i = 0; i < p; i++)
        D[i] = new double* [q];
    for (int i = 0; i < p; i++)
        for (int j = 0; j < q; j++)
            D[i][j] = new double[m];

    C = new double* [p];
    for (int i = 0; i < p; i++)
        C[i] = new double[q];

    print_two_dimensional_matrixs();
    calculate_matrix_C();
}

void SIMDMatrixes::print_two_dimensional_matrixs()
{
    cout << endl << "Матрица A" << endl << endl;
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < m; j++)
            cout << setw(7) << A[i][j] << "\t";
        cout << endl;
    }

    cout << endl << "Матрица B" << endl << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < q; j++)
            cout << setw(7) << B[i][j] << "\t";
        cout << endl;
    }

    cout << endl << "Матрица E" << endl << endl;
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < m; j++)
            cout << setw(7) << E[i][j] << "\t";
        cout << endl;
    }

    cout << endl << "Матрица G" << endl << endl;
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++)
            cout << setw(7) << G[i][j] << "\t";
        cout << endl;
    }
}

void SIMDMatrixes::calculate_matrix_C()
{
    Fijk(A, B, E, p, m, q);
    Dijk(A, B, p, m, q);
  
    Cij(F, G, D, p, q, m);
    print_auxiliary_matrix();
    print_final_matrix();
}

void SIMDMatrixes::print_auxiliary_matrix()
{
    cout << endl << "Вспомогательная матрица Fij" << endl << endl;

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            for (int k = 0; k < m; k++) {
                cout << setw(7) << F[i][j][k];
            }
            cout << endl;
        }
        cout << endl;
    }

    cout << endl << "Вспомогательная матрица Dij" << endl << endl;

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            for (int k = 0; k < m; k++) {
                cout << setw(7) << F[i][j][k];
            }
            cout << endl;
        }
        cout << endl;
    }
}

void SIMDMatrixes::print_final_matrix()
{
    cout << "Матрица C:";
    for (int i = 0; i < p; i++) {
        cout << endl;
        for (int j = 0; j < q; j++) {
            cout << setw(12) << C[i][j] << " ";
        }
    }
}

void SIMDMatrixes::Cij(double*** F, double** G, double*** D, int i, int j, int k) {
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            C[i][j] = singleConjunction(i,j) * (3. * G[i][j] - 2.) * G[i][j] + (singleDisjunction(i, j) + (4. * multiplication(i, j) - 3. * singleDisjunction(i, j)) * G[i][j]) * (1. - G[i][j]);
        }
    }
}

void SIMDMatrixes::Fijk(double** A, double** B, double** E, int i, int j, int k) {
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < m; k++)
            {
                F[i][j][k] = (contraction(A[i][k],B[k][j])) * (2. * E[0][k] - 1.) * E[0][k] + uncontraction(A, B, i, k, j) * (1. + (4. * contraction(A[i][k], B[k][j]) - 2.) * E[0][k]) * (1. - E[0][k]);
            }
        }
    }
 }

void SIMDMatrixes::Dijk(double** A, double** B, int i, int j, int k) {
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < m; k++)
            {
                D[i][j][k] = conjunction_self(A[i][k], B[k][j]);
            }
        }
    }
}

double SIMDMatrixes::singleConjunction(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        result *= F[i][j][k];
    }
    return result;
}

double SIMDMatrixes::singleDisjunction(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        result *= 1 - D[i][j][k];
    }
    
    return 1 - result;
}

double SIMDMatrixes::multiplication(int i, int j) {
    double res, half;
    half = (singleConjunction(i, j) + singleDisjunction(i, j) - 1.);
    if (half > 0.) res = half;
    else res = 0.;
 
    return res;
}

double SIMDMatrixes::contraction(int first,int second) {
    int delta = 1;
    while (min((1 - first), delta) > second)
    {
        delta -= 1;

    }
    return delta;

}

double SIMDMatrixes::uncontraction(double** A, double** B, int i, int k, int j) {
    double res;
    if ((1 - B[k][j]) > A[i][k]) res = (1 - B[k][j]);
    else res = A[i][k];
  
    return res;
}

double SIMDMatrixes::conjunction_self(int first,int second) {
    return min(first, second);
}