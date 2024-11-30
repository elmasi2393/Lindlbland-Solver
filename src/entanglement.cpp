#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <complex>
#include <Qobj_eigen.cpp>
#include <chrono>
// #include <utils.cpp>    // Sobrecarga de operadores necesarios
// #include <Qobj.cpp>     // Definición de matrices y funciones de matrices

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>

#define hbar 1.0

using namespace std;
using namespace Eigen;


// ---------- Operadores evolucion temporal ------------------

class Hamiltonian : public QobjEvo {
    public:
        Hamiltonian(const double wq1, const double wq2) : wq1(wq1), wq2(wq2), operator_a {tensor_product(sigma_z, sigma_0)}, operator_b {tensor_product(sigma_0, sigma_z)} {}
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return wq1/2 * operator_a + wq2/2 * operator_b;
        }
    private:
        double wq1, wq2;
        QuantumMatrix operator_a, operator_b;
};

class C_phi : public QobjEvo {
    public:
        C_phi(const double gamma_min_1, const double gamma_min_2) : gamma_min_1(gamma_min_1), gamma_min_2(gamma_min_2), operator_a(tensor_product(sigma_z, sigma_0)) , operator_b(tensor_product(sigma_0, sigma_z)) {}
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return sqrt(gamma_min_1) * operator_a + sqrt(gamma_min_2) * operator_b;
        }
    private:
        double gamma_min_1, gamma_min_2;
        QuantumMatrix operator_a, operator_b;
};

// class C_phi_global : public QobjEvo {
//     public:
//         C_phi_global(const double gamma_min) : gamma_min(gamma_min), operator_a(tensor_product(sigma_z, sigma_0)) , operator_b(tensor_product(sigma_0, sigma_z)) {}
//         QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
//             return sqrt(gamma_min) * (operator_a + operator_b);
//         }
//     private:
//         double gamma_min;
//         QuantumMatrix operator_a, operator_b;
// };

// class C_minnus : public QobjEvo {
//     public:
//         C_minnus(const double gamma_min_1, const double gamma_min_2) : gamma_min_1(gamma_min_1), gamma_min_2(gamma_min_2), operator_a(tensor_product(sigma_m, sigma_0)) , operator_b(tensor_product(sigma_0, sigma_m)) {}
//         QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
//             return sqrt(gamma_min_1) * operator_a + sqrt(gamma_min_2) * operator_b;
//         }
//     private:
//         double gamma_min_1, gamma_min_2;
//         QuantumMatrix operator_a, operator_b;
// };


int main(){
    double wq1 = 1.0;
    double wq2 = 1.0;
    double gamma_phi_1 = 0.1;
    double gamma_phi_2 = 0.1;

    // QuantumVector phi_plus = bell_state("00");

    vector<double> diag = {1.0/3.0, 1.0/6.0, 1.0/6.0, 1.0/3.0};
    QuantumValue w = 1.0/3.0;
    QuantumValue z = 0.0;
    QuantumMatrix density_m = standar_death_form(diag, w, z);


    string path2save = "../results/results_standar_death_form.csv";
    
    QobjEvo * C_p = new C_phi(gamma_phi_1, gamma_phi_2);
    QobjEvo * H = new Hamiltonian(wq1, wq2);

    vector<QobjEvo*> C_ops = {C_p};

    // // Definición de los tiempos
    double h = 0.001;
    vector<double> time_limits = {0, 2};

    cout << "Resolviendo la ecuación de Lindblad" << endl;
    cout << "Tiempo inicial: " << time_limits[0] << endl;
    cout << "Tiempo final: " << time_limits[1] << endl;
    cout << "Paso de tiempo: " << h << endl;
    cout << "Cantidad de C_ops: " << C_ops.size() << endl;
    cout << "Tamaño de la matriz de densidad: " << density_m.rows() << "x" << density_m.cols() << endl;

    
    // // Resuelvo la ecuación de Lindblad
    auto start = chrono::high_resolution_clock::now();
    vector<QuantumMatrix> states = mesolve(H, time_limits, h, density_m, C_ops);
    auto end = chrono::high_resolution_clock::now();

    cout << "Resolución de la ecuación de Lindblad completada" << endl;
    cout << "Tiempo de ejecución: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "s" << endl;

    // // Guardo los resultados
    cout << "Guardando los resultados en el archivo: " << path2save << endl;
    save_states(states, h, time_limits, path2save);
    cout << "Resultados guardados" << endl;
    cout << "Proceso completado" << endl;

    return 0;
}