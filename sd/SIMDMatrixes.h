#pragma once

using namespace std;

class SIMDMatrixes {

private:

    int m;
    int p;
    int q;

    int process_element;

    int r;

    double** E;

    double** A;
    double** B;
    double** G;


    double*** F;
    double*** D;

    double** C;


    int T1;
    int Tn;
    double Ky;
    double Eff;
    double Diff;
    int Lavg;
    
public:
    
        double get_random_elements();
        void filling_matrixs();
        void print_two_dimensional_matrixs();
        void calculate_matrix_C();
        void print_auxiliary_matrix();
        void set_input_elements_of_massiv();



        double cij(int i, int j);

        void Fijk(double** A, double** B, double** E, int i, int j, int k);
        void Dijk(double** A, double** B, int i, int j, int k);


        double a_and_b(double** A, double** B, int i, int k, int j);

        double a_to_b(int first, int second);

        double b_to_a(double** A, double** B, int i, int k, int j);

        double f_func(int i, int j);

        double d_func(int i, int j);

        double f_and_d(int i, int j);

};
