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
        Hamiltonian(const double wq) : wq(wq) {}
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return wq/2 * sigma_z;
        }
    private:
        double wq;
};

class C_minnus : public QobjEvo {
    public:
        C_minnus(const double gamma_min): gamma_min(gamma_min) {}
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return sqrt(gamma_min) * sigma_m;
        }
    private:
        double gamma_min;
};

class C_plus : public QobjEvo {
    public:
        C_plus(const double gamma_plus): gamma_plus(gamma_plus) {}
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return sqrt(gamma_plus) * sigma_p;
        }
    private:
        double gamma_plus;
};

class C_phi : public QobjEvo {
    public:
        C_phi(const double gamma_phi): gamma_phi(gamma_phi) {};
        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const override {
            return sqrt(gamma_phi/2) * sigma_z;
        }
    private:
        double gamma_phi;
};


void  run_simulation(const vector<double> system_parameters, const vector<double> time_limits, const double step, const QuantumMatrix &p0, const string path2save){
    const double wq = system_parameters[0];
    const double gamma_min = system_parameters[1];
    const double gamma_plus = system_parameters[2];
    const double gamma_phi = system_parameters[3];

    QobjEvo * H = new Hamiltonian(wq);
    QobjEvo * C_m = new C_minnus(gamma_min);
    QobjEvo * C_p = new C_plus(gamma_plus);
    QobjEvo * C_ph = new C_phi(gamma_phi);
    vector<QobjEvo*> C_ops = {C_m, C_p, C_ph};

    cout << "Parámetros del sistema:" << endl;
    cout << "{wq, gamma_min} = {" << wq << ", " << gamma_min << ", " << gamma_plus << ", " << gamma_phi << "}" << endl;
    cout << "Tiempo inicial: " << time_limits[0] << endl;
    cout << "Tiempo final: " << time_limits[1] << endl;
    cout << "Paso de tiempo: " << step << endl;
    cout << "Cantidad de C_ops: " << C_ops.size() << endl;
    cout << "Tamaño de la matriz de densidad: " << p0.rows() << "x" << p0.cols() << endl;
    cout << "Estado inicial: " << endl;
    cout << p0 << endl;

    cout << "Resolviendo la ecuación de Lindblad" << endl;
    auto start = chrono::high_resolution_clock::now();
    vector<QuantumMatrix> states = mesolve(H, time_limits, step, p0, C_ops);
    auto end = chrono::high_resolution_clock::now();
    cout << "Tiempo de ejecución: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "s" << endl;
    cout << "Guardando los estados en el archivo: " << path2save << endl;
    save_states(states, step, time_limits, path2save);
    cout << "Guardado completado" << endl;
    cout << "Programa finalizado" << endl;

    return;
}

int main(){
    double wq = 1.0;
    double gamma_min = 0.1;
    double gamma_plus = 0.1;
    double gamma_phi = 0.1;

    vector<double> parameters = {wq, gamma_min, gamma_plus, gamma_phi};
    string path2save = "../results/results_psi_plus_relax_exc_desf.csv";

    // vector<QuantumVector> bell_states = {bell_state("00"), bell_state("01"), bell_state("10"), bell_state("11")};
    
    // QuantumMatrix p0 = ket2dm(x_plus_state);
    QuantumVector p = (QuantumVector(2) << sqrt(3.0/4.0), sqrt(1.0/4.0)).finished();
    QuantumMatrix p0 = ket2dm(p);


    vector<double> time_limits = {0, 30};
    double h = 0.001;

    run_simulation(parameters, time_limits, h, p0, path2save);
    
    
    return 0;
}