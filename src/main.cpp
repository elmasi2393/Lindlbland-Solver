/* ------------------ TRABAJO FINAL ------------------

Resolución de la ecuación de Lindbland:

dp/dt = -i[H, p] + 1/2 sum_n [2 Cn p Cn^d - Cn^d Cn p - p Cn^d Cn]

donde H es el Hamiltoniano, p es la matriz densidad y Cn son los operadores de colapso.

Para resolver la ecuación de Lindbland, se puede utilizar el método de Runge-Kutta de cuarto orden. Para ello, se debe definir un functor que contenga la derivada de la matriz densidad, y un solver que resuelva la ecuación diferencial.
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <complex>

#include <utils.cpp>    // Sobrecarga de operadores necesarios
#include <Qobj.cpp>     // Definición de matrices y funciones de matrices

using namespace std;


// Defino los operadores de colapso y hamiltoniano
// class C_minnus : public QobjEvo {
//     public:
//         C_minnus(const double gamma_min) : gamma_min(gamma_min) {}
//         quantum_object operator()(const double t, const quantum_object &variables) const override {
//             return sqrt(gamma_min) * sigma_m;
//         }
//     private:
//         double gamma_min;
// };

// class C_plus : public QobjEvo {
//     public:
//         C_plus(const double gamma_plus) : gamma_plus(gamma_plus) {}
//         quantum_object operator()(const double t, const quantum_object &variables) const override {
//             return sqrt(gamma_plus) * sigma_p;
//         }
//     private:
//         double gamma_plus;
// };

// class C_phi : public QobjEvo {
//     public:
//         C_phi(const double gamma_phi): gamma_phi(gamma_phi) {}
//         quantum_object operator()(const double t, const quantum_object &variables) const override{
//             return sqrt(gamma_phi/2) * sigma_z;
//         }
//     private:
//         double gamma_phi;
// };

// class Hamiltonian : public QobjEvo {
//     public:
//         Hamiltonian(const double wq) : wq(wq) {}
//         quantum_object operator()(const double t, const quantum_object &variables) const override {
//             return wq/2 * sigma_z;
//         }
//     private:
//         double wq;
// };

class Hamiltonian_2qubit : public QobjEvo {
    public:
        Hamiltonian_2qubit(const double wq) : wq(wq) {}
        quantum_object operator()(const double t, const quantum_object &variables) const override {
            return wq/2 * tensor(sigma_z, sigma_z);
        }
    private:
        double wq;
};

class C_phi_2qubit : public QobjEvo{
    public:
        C_phi_2qubit(const double gamma_phi) : gamma_phi(gamma_phi) {}
        quantum_object operator()(const double t, const quantum_object &variables) const override {
            return sqrt(gamma_phi/2) * tensor(sigma_z, sigma_z);
        }
    private:
        double gamma_phi;
};

int main(){
    double wq = 1.0;
    double gamma_min = 0.1;
    double gamma_plus = 0.1;
    double gamma_phi = 0.1;
    quantum_object x_state = basis(2, 0) + basis(2, 1);

    string path2save = "../results/results_x_state.csv";

    // QobjEvo * C_m = new C_minnus(gamma_min);
    // QobjEvo * C_p = new C_plus(gamma_plus);
    QobjEvo * C_ph = new C_phi(gamma_phi);

    QobjEvo * H = new Hamiltonian(wq);

    vector<QobjEvo*> C_ops = {C_ph};

    // Definición de los tiempos
    double h = 0.0001;
    vector<double> time_limits = {0, 30};
    vector<double> time_axis((int)time_limits[1]/h + 1, 0);

    for (int i = 0; i < time_axis.size(); i++)
        time_axis[i] = i * h;

    cout << "Time axis: " << time_axis.size() << endl;
    cout << "Time limits: " << time_axis[0] << "," << time_axis[time_axis.size() - 1] << endl;

    // Resuelvo la ecuación de Lindblad
    // vector<quantum_object> states = mesolve(H, time_axis, normalize(x_state), C_ops);

    cout << tensor(sigma_z, sigma_x) << endl;
    // cout << "Resultados: " << endl;
    // cout << states;

    // Guardo los resultados
    // save_states(states, h, time_limits, path2save);

    return 0;
}