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
    int T1;
    int Tn;
    double Ky;
    double Eff;
    double Diff;
    int Lavg;

    //**результатирующая матрица**//
    double** C;
    
public:
    
        double get_random_elements();

        void set_input_elements_of_massiv();
        void filling_matrixs();
        void calculate_matrix_C();
        void print_two_dimensional_matrixs();
        void print_auxiliary_matrix();
        void print_final_matrix();


        //**Нахождение матриц F,D,C**//
        void Fijk(double** A, double** B, double** E, int i, int j, int k);
        void Dijk(double** A, double** B, int i, int j, int k);
        
        void Cij(double*** F, double** G, double*** D, int i, int j, int k);

        //**Вычисление дополнитеьных матриц Fij и Dij
        double contraction(int first, int second);
        double conjunction_self(int first, int second);

        double uncontraction(double** A, double** B, int i, int k, int j);
        double singleConjunction(int i, int j);
        double singleDisjunction(int i, int j);
        double multiplication(int i, int j);

};
