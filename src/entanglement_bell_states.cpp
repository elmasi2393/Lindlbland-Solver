#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <complex>
#include <Qobj_eigen.cpp>
#include <chrono>
#include <thread>

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


void  run_simulation(const vector<double> system_parameters, const vector<double> time_limits, const double step, const QuantumMatrix &p0, const string path2save){
    const double wq1 = system_parameters[0];
    const double wq2 = system_parameters[1];
    const double gamma_phi_1 = system_parameters[2];
    const double gamma_phi_2 = system_parameters[3];

    QobjEvo * H = new Hamiltonian(wq1, wq2);
    QobjEvo * C_p = new C_phi(gamma_phi_1, gamma_phi_2);
    vector<QobjEvo*> C_ops = {C_p};

    cout << "Parámetros del sistema:" << endl;
    cout << "{wq1, wq2, gamma_phi_1, gamma_phi_2} = {" << wq1 << ", " << wq2 << ", " << gamma_phi_1 << ", " << gamma_phi_2 << "}" << endl;
    cout << "Tiempo inicial: " << time_limits[0] << endl;
    cout << "Tiempo final: " << time_limits[1] << endl;
    cout << "Paso de tiempo: " << step << endl;
    cout << "Cantidad de C_ops: " << C_ops.size() << endl;
    cout << "Tamaño de la matriz de densidad: " << p0.rows() << "x" << p0.cols() << endl;
    cout << "Estado inicial: " << endl;
    cout << p0 << endl;

    cout << "Resolviendo la ecuación de Lindblad" << endl;
    auto start = chrono::high_resolution_clock::now();
    vector<QuantumMatrix> states = mesolve(H, time_limits, step, ket2dm(p0), C_ops);
    auto end = chrono::high_resolution_clock::now();
    cout << "Tiempo de ejecución: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "s" << endl;
    cout << "Guardando los estados en el archivo: " << path2save << endl;
    save_states(states, step, time_limits, path2save);
    cout << "Guardado completado" << endl;
    cout << "Programa finalizado" << endl;

    return;
}

int main(){
    // Definicion de variables del sistema y operadores de colapso
    double wq1 = 1.0;
    double wq2 = 1.0;
    double gamma_phi_1 = 0.1;
    double gamma_phi_2 = 0.1;

    // Acomodo los parámetros y los estados de Bell en vectores para pasarle a la simulación
    vector<double> parameters = {wq1, wq2, gamma_phi_1, gamma_phi_2}; 
    vector<QuantumVector> bell_states = {bell_state("00"), bell_state("01"), bell_state("10"), bell_state("11")};

    // defino escala temporal
    vector<double> time_limits = {0, 10};
    double h = 0.001;

    vector<thread> threads; // hilos para ejecución en paralelo

    for(int i = 0; i < bell_states.size(); i++){    // crea los hilos y la ruta donde guardar los resulatdos
        string path2save = "../results/results_bell_state_" + to_string(i/2) + to_string(i%2) + ".csv";
        threads.push_back(thread(run_simulation, parameters, time_limits, h, bell_states[i], path2save));
    }

    // Espero a que terminen
    for(auto &t : threads)
        t.join();
    
    return 0;
}