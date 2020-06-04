#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>

#include "SIMDMatrixes.h"

using namespace std;


int main() {
    
        
        srand(time(nullptr));
        setlocale(LC_ALL, "Russian");

        SIMDMatrixes lab2;

        lab2.set_input_elements_of_massiv();
        lab2.filling_matrixs();
        

   
}

