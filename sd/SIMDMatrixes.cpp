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
        else {      
            filling_matrixs();
            break;
        }
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
    Fijk();
    Dijk();
    Cij();
    print_auxiliary_matrix();
    print_final_matrix();

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

void SIMDMatrixes::print_final_matrix()
{
    cout << "������� C:";
    for (int i = 0; i < p; i++) {
        cout << endl;
        for (int j = 0; j < q; j++) {
            cout << setw(12) << C[i][j] << " ";
        }
    }
    print_elements();
}

void SIMDMatrixes::print_elements()
{
    rang = m * p * q;
    T1 = numberOfMultiplications * multiplicationTime + numberOfAdditions * additionTime +
        numberOfSubtractions * subtractionTime + numberOfComparisons * comparisonTime;
    Ky = T1 / Tn;
    e = Ky / process_element;
    Lavg = ceil(averageTime / rang);
    diff = (Tn / Lavg);
    cout << "\n\n���� = " << rang << "\nT1 = " << T1 << "\nTn = " << Tn << "\nKy = " << Ky << "\ne = " << e << "\nLavg = " << Lavg << "\nDiff = " << diff << "\nLsum = " << averageTime;
}

void SIMDMatrixes::Cij() {
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            C[i][j] = singleConjunction(i,j) * (3 * G[i][j] - 2) * G[i][j] + (singleDisjunction(i, j) + (4 * multiplication(i, j) - 3 * singleDisjunction(i, j)) * G[i][j]) * (1 - G[i][j]);
            numberOfMultiplications += 8;
            numberOfAdditions += 3;
            numberOfSubtractions += 2;
        }
    }
    Tn += ceil((p * q) / process_element) * (6 * multiplicationTime + 2 * additionTime + 2 * subtractionTime + 3 * ((m - 1) * multiplicationTime + (m + 1) * subtractionTime) +
        2 * ((m - 1) *
            multiplicationTime) +
        subtractionTime + additionTime +
        2 * comparisonTime);
    averageTime += p * q * (6 * multiplicationTime + 2 * additionTime + 2 * subtractionTime + 3 * ((m - 1) * multiplicationTime + (m + 1) * subtractionTime) + 2 * ((m - 1) * multiplicationTime) + subtractionTime + additionTime + 2 * comparisonTime);
}

void SIMDMatrixes::Fijk() {
    for (int i = 0; i < p; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < m; k++)
            {
                F[i][j][k] = (contraction(A[i][k], B[k][j])) * (2 * E[0][k] - 1) * E[0][k] + uncontraction(i, k, j) * (1 + (4 * contraction(A[i][k], B[k][j]) - 2) * E[0][k]) * (1 - E[0][k]);
                numberOfMultiplications += 7;
                numberOfAdditions += 2;
                numberOfSubtractions += 3;
            }
        }
    }
    Tn += ceil((p * q * m) / process_element) * (7 * multiplicationTime + 2 * additionTime + 3 * subtractionTime +
        3 * (2 * comparisonTime + subtractionTime));
    averageTime += p * q * m * (7 * multiplicationTime + 2 * additionTime + 3 * subtractionTime + 3 * (2 * comparisonTime + subtractionTime));
}

void SIMDMatrixes::Dijk() {
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
    Tn += ceil((p * q * m) / process_element) * 2 * comparisonTime;
    averageTime += p * q * m * 2 * comparisonTime;
}

double SIMDMatrixes::singleConjunction(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        result *= F[i][j][k];
        numberOfMultiplications++;
    }
    numberOfMultiplications--;
    return result;
}

double SIMDMatrixes::singleDisjunction(int i, int j) {
    double result = 1;
    for (int k = 0; k < m; k++) {
        numberOfMultiplications++;
        numberOfSubtractions++;
        result *= 1 - D[i][j][k];
    }
    numberOfMultiplications--;
    numberOfSubtractions++;
    return 1 - result;
}

double SIMDMatrixes::multiplication(int i, int j) {
    double res, half;
    numberOfAdditions++;
    numberOfSubtractions++;
    numberOfComparisons += 2;
    half = (singleConjunction(i, j) + singleDisjunction(i, j) - 1,0);
    if (half > 0) res = half;
    else res = 0;
 
    return res;
}

double SIMDMatrixes::contraction(int first,int second) {
    int delta = 1;
    while (min((1 - first), delta) > second)
    {
        delta -= 1;
        numberOfComparisons += 2;
        numberOfSubtractions += 2;

    }
    return delta;

}

double SIMDMatrixes::uncontraction(int i, int k, int j) {
    double res;
    if ((1 - B[k][j]) > A[i][k]) {
        res = (1 - B[k][j]);
        numberOfComparisons += 2;
        numberOfSubtractions += 2;
    }
    else res = A[i][k];
  
    return res;
}

double SIMDMatrixes::conjunction_self(int first,int second) {
    numberOfComparisons += 2;
    return min(first, second);
}