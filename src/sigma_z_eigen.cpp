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
#define hbar 1.0

using namespace std;
using namespace Eigen;


// ---------- Operadores evolucion temporal ------------------

class Hamiltonian : public QobjEvo {
    public:
        Hamiltonian(const double wq) : wq(wq) {}
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return wq/2 * sigma_z;
        }
    private:
        double wq;
};

class C_minnus : public QobjEvo {
    public:
        C_minnus(const double gamma_min) : gamma_min(gamma_min) {}
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return sqrt(gamma_min) * sigma_m;
        }
    private:
        double gamma_min;
};



int main(){
    double wq = 1.0;
    double gamma_min = 0.1;
    double gamma_plus = 0.1;
    double gamma_phi = 0.1;

    // quantum_object exc_state = tensor(basis(2, 0), basis(2, 0));
    QuantumVector exc = basis(2, 0);
    QuantumMatrix density_m = exc * exc.adjoint();


    string path2save = "../results/results_eigen_sigmaz.csv";
    
    QobjEvo * C_m = new C_minnus(gamma_min);
    QobjEvo * H = new Hamiltonian(wq);

    vector<QobjEvo*> C_ops = {C_m};

    // // Definición de los tiempos
    double h = 0.001;
    vector<double> time_limits = {0, 30};

    
    // // Resuelvo la ecuación de Lindblad
    auto start = chrono::high_resolution_clock::now();
    vector<QuantumMatrix> states = mesolve(H, time_limits, h, density_m, C_ops);
    auto end = chrono::high_resolution_clock::now();

    cout << "Tiempo de ejecución: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "s" << endl;


    // // Guardo los resultados
    save_states(states, h, time_limits, path2save);

    return 0;
}