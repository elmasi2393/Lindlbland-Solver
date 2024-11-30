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
#define hbar 1.0

using namespace std;


class Hamiltonian : public QobjEvo {
    public:
        Hamiltonian(const double wq_1, const double wq_2) : wq_1(wq_1), wq_2(wq_2) {}
        quantum_object operator()(const double t, const quantum_object &variables) const override {
            return wq_1/2 * tensor(sigma_z, I) + wq_2/2 * tensor(I, sigma_z);
        }
    private:
        double wq_1, wq_2;
};

class C_minnus : public QobjEvo {
    public:
        C_minnus(const double gamma_min_1, const double gamma_min_2) : gamma_min_1(gamma_min_1), gamma_min_2(gamma_min_2) {}
        quantum_object operator()(const double t, const quantum_object &variables) const override {
            return sqrt(gamma_min_1) * tensor(sigma_m, I);
        }
    private:
        double gamma_min_1, gamma_min_2;
};

// class C_phi_2qubit : public QobjEvo{
//     public:
//         C_phi_2qubit(const double gamma_phi) : gamma_phi(gamma_phi) {}
//         quantum_object operator()(const double t, const quantum_object &variables) const override {
//             return sqrt(gamma_phi/2) * tensor(sigma_z, sigma_z);
//         }
//     private:
//         double gamma_phi;
// };

int main(){
    double wq = 1.0;
    double gamma_min = 0.1;
    double gamma_plus = 0.1;
    double gamma_phi = 0.1;
    quantum_object exc_state = tensor(basis(2, 0), basis(2, 0));
    

    string path2save = "../results/results_2q_sigmaz_1.csv";
    
    QobjEvo * C_m = new C_minnus(gamma_min, gamma_min);
    QobjEvo * H = new Hamiltonian(wq, wq);

    vector<QobjEvo*> C_ops = {C_m};

    // Definición de los tiempos
    double h = 0.001;
    vector<double> time_limits = {0, 30};
    vector<double> time_axis((int)time_limits[1]/h + 1, 0);

    for (int i = 0; i < time_axis.size(); i++)
        time_axis[i] = i * h;

    cout << "Time axis: " << time_axis.size() << endl;
    cout << "Time limits: " << time_axis[0] << "," << time_axis[time_axis.size() - 1] << endl;

    // Resuelvo la ecuación de Lindblad
    vector<quantum_object> states = mesolve(H, time_axis, normalize(exc_state), C_ops);

    cout << exc_state;
    // cout << tensor(I, I) * exc_state;
    // cout << tensor(sigma_z, I) * exc_state;
    // Guardo los resultados
    save_states(states, h, time_limits, path2save);

    return 0;
}