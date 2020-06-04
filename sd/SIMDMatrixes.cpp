#include "SIMDMatrixes.h"
#include <iostream>
#include <iomanip>
#include <array>



using namespace std;

void SIMDMatrixes::set_input_elements_of_massiv()
{
    while (true)
    {
        cout << "������� m: " << endl;
        cin >> m;
        cout << "������� p: " << endl;
        cin >> p;
        cout << "������� q: " << endl;
        cin >> q;
        cout << "������� ���������� ����-�� ���������: " << endl;
        cin >> process_element;
        if (m == 0 || p == 0 || q == 0 || process_element == 0) {
            cout << "������!��������� ���� ������\n";
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

    print_two_dimensional_matrixs();
    calculate_matrix_C();
}

void SIMDMatrixes::print_two_dimensional_matrixs()
{
    cout << endl << "������� A" << endl << endl;
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < m; j++)
            cout << setw(7) << A[i][j] << "\t";
        cout << endl;
    }

    cout << endl << "������� B" << endl << endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < q; j++)
            cout << setw(7) << B[i][j] << "\t";
        cout << endl;
    }

    cout << endl << "������� E" << endl << endl;
    for (int i = 0; i < 1; i++) {
        for (int j = 0; j < m; j++)
            cout << setw(7) << E[i][j] << "\t";
        cout << endl;
    }

    cout << endl << "������� G" << endl << endl;
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
    print_auxiliary_matrix();

}

void SIMDMatrixes::print_auxiliary_matrix()
{
    cout << endl << "��������������� ������� Fij" << endl << endl;

    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            for (int k = 0; k < m; k++) {
                cout << setw(7) << F[i][j][k];
            }
            cout << endl;
        }
        cout << endl;
    }

    cout << endl << "��������������� ������� Dij" << endl << endl;

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



//deafult task c(ij)
double SIMDMatrixes::cij(int i, int j) {
    double c;
   
    c = f_func(i, j) * (3. * G[i][j] - 2.) * G[i][j] + (d_func(i, j) + (4. * f_and_d(i, j) - 3. * d_func(i, j)) * G[i][j]) * (1. - G[i][j]);

    return c;
}

//deafult task f(ijk)
void SIMDMatrixes::Fijk(double** A, double** B, double** E, int i, int j, int k) {
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < m; k++)
            {

                F[i][j][k] = (a_to_b(A[i][k],B[k][j])) * (2. * E[0][k] - 1.) * E[0][k] + b_to_a(A, B, i, k, j) * (1. + (4. * a_to_b(A[i][k], B[k][j]) - 2.) * E[0][k]) * (1. - E[0][k]);
            }
        }
    }
 }

//deafult task d(ijk)
void SIMDMatrixes::Dijk(double** A, double** B, int i, int j, int k) {
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < m; k++)
            {
                D[i][j][k] = a_and_b(A, B, i, k, j);
            }
        }
    }
}

//F func
double SIMDMatrixes::f_func(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        result *= F[i][j][k];
    }
    return result;
}


//D func
double SIMDMatrixes::d_func(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        result *= 1 - D[i][j][k];
    }
    
    return 1 - result;
}




//F and D
double SIMDMatrixes::f_and_d(int i, int j) {
    double res, half;
    half = (f_func(i, j) + d_func(i, j) - 1.);
    if (half > 0.) res = half;
    else res = 0.;
 
    return res;
}


//A to B
double SIMDMatrixes::a_to_b(int first,int second) {
    int delta = 1;
    while (min((1 - first), delta) > second)
    {
        delta -= 1;

    }
    return delta;

}


//B to A
double SIMDMatrixes::b_to_a(double** A, double** B, int i, int k, int j) {
    double res;
    if ((1 - B[k][j]) > A[i][k]) res = (1 - B[k][j]);
    else res = A[i][k];
  
    return res;
}


//A and B
double SIMDMatrixes::a_and_b(double** A, double** B, int i, int k, int j) {
    double res;
    if (A[i][k] > B[k][j]) res = A[i][k];
    else res = B[k][j];
    return res;
}