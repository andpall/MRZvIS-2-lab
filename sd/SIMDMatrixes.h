#pragma once

using namespace std;

class SIMDMatrixes {

private:

    //**размеры матриц**//
    int m;
    int p;
    int q;

    //**элемент который влияет на работу**//
    int process_element;

    //**ранг**//
    int r;

    //**двухмерные матрицы**//
    double** E;
    double** A;
    double** B;
    double** G;

    //**трехмерные матрицы**//
    double*** F;
    double*** D;


    //**элементы для подсчёта**//
    int multiplicationTime = 1;
    int additionTime = 1;
    int subtractionTime = 1;
    int comparisonTime = 1;

    int rang;
    double Ky;
    double e;
    double Lavg;
    double Tn;
    double T1;
    double diff;

    int numberOfMultiplications = 0;
    int numberOfAdditions = 0;
    int numberOfSubtractions = 0;
    int numberOfComparisons = 0;
    int averageTime = 0;

    //**результатирующая матрица**//
    double** C;
    
public:
    
        double get_random_elements();

        void set_input_elements_of_massiv();
        void filling_matrixs();
        void calculate_matrix_C();
        void print_two_dimensional_matrixs();
        void print_final_matrix();

        void print_elements();

        //**Нахождение матриц F,D,C**//
        void Fijk();
        void Dijk();
        void Cij();

        //**Вычисление
        double contraction(int first, int second);
        double conjunction_self(int first, int second);
        double uncontraction(int i, int k, int j);
        double singleConjunction(int i, int j);
        double singleDisjunction(int i, int j);
        double multiplication(int i, int j);

        void print_auxiliary_matrix();

};
