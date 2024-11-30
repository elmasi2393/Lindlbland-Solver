#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Definimos un tipo para matrices complejas en Eigen
using QuantumMatrix = MatrixXcd;  // Matriz de números complejos de tamaño dinámico
using QuantumVector = VectorXcd;  // Vector de números complejos de tamaño dinámico


int main() {
    // Crear una matriz cuántica compleja de ejemplo (4x1)
    QuantumMatrix matrix(2, 2);
    matrix << 1/sqrt(2), 0,
              0, 1/sqrt(2);

    QuantumVector state(2);
    state << 1/sqrt(2), 1/sqrt(2);

    cout << matrix * state << endl;
    cout << state.rows() << endl;
    cout << state.cols() << endl;

    cout << matrix.adjoint() << endl;
    return 0;

}
