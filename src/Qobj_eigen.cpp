
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <complex>

// #include <utils.cpp>    // Sobrecarga de operadores necesarios
// #include <Qobj.cpp>     // Definición de matrices y funciones de matrices

#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#define hbar 1.0

using namespace std;
using namespace Eigen;

// Definimos un tipo para matrices complejas en Eigen
using QuantumMatrix = MatrixXcd;  // Matriz de números complejos de tamaño dinámico
using QuantumVector = VectorXcd;  // Vector de números complejos de tamaño dinámico
using QuantumValue = complex<double>;


// -------------------------- Definicion matrices de Pauli ---------------------------------------
static const QuantumMatrix sigma_x = (QuantumMatrix(2, 2) << 0., 1., 1., 0.).finished();
static const QuantumMatrix sigma_y = (QuantumMatrix(2, 2) << 0., QuantumValue(0, -1), QuantumValue(0, 1), 0.).finished();
static const QuantumMatrix sigma_z = (QuantumMatrix(2, 2) << 1., 0., 0., -1.).finished();
static const QuantumMatrix sigma_m = (QuantumMatrix(2, 2) << 0., 0., 1., 0.).finished();
static const QuantumMatrix sigma_p = (QuantumMatrix(2, 2) << 0., 1., 0., 0.).finished();
static const QuantumMatrix sigma_0 = (QuantumMatrix(2, 2) << 1., 0., 0., 1.).finished();

// -------------------------- Definicion de estados de un TLS --------------------------------------
static const QuantumVector z_up_state = (QuantumVector(2) << 1., 0.).finished();
static const QuantumVector z_down_state = (QuantumVector(2) << 0., 1.).finished();
static const QuantumVector x_plus_state = (QuantumVector(2) << 1./sqrt(2), 1./sqrt(2)).finished();
static const QuantumVector x_minus_state = (QuantumVector(2) << 1./sqrt(2), -1./sqrt(2)).finished();
static const QuantumVector y_plus_state  = (QuantumVector(2) << 1./sqrt(2), QuantumValue(0, 1./sqrt(2))).finished();
static const QuantumVector y_minus_state = (QuantumVector(2) << 1./sqrt(2), QuantumValue(0, -1./sqrt(2))).finished();


// --------------------- Operador de evolución temporal -------------------------------------------
class QobjEvo {
    public:
        virtual QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const = 0;
};


// -------------------- Operaciones sobre estados/operadores --------------------------------------
// Conmutador entre operadores
QuantumMatrix conmutator(const QuantumMatrix &q1, const QuantumMatrix &q2){
    return q1 * q2 - q2 * q1;
}

// Función para el producto tensorial
QuantumMatrix tensor_product(const QuantumMatrix& A, const QuantumMatrix& B) {
    // int rows_A = A.rows();
    // int cols_A = A.cols();
    // int rows_B = B.rows();
    // int cols_B = B.cols();

    // // Crear una matriz de tamaño (rows_A * rows_B) x (cols_A * cols_B)
    // QuantumMatrix result(rows_A * rows_B, cols_A * cols_B);

    // // Rellenar la matriz resultante
    // for (int i = 0; i < rows_A; ++i) {
    //     for (int j = 0; j < cols_A; ++j) {
    //         // El valor A(i, j) se multiplica por toda la matriz B
    //         result.block(i * rows_B, j * cols_B, rows_B, cols_B) = A(i, j) * B;
    //     }
    // }
    // return result;
    return kroneckerProduct(A, B);
}


// ---------------------- Ecuación de Lindblad para el Solver -----------------------------------------------
class LindbladEvo : public QobjEvo {
    public:
        LindbladEvo(QobjEvo * H, const vector<QobjEvo*> &C_ops) : H(H), C_ops(C_ops) {}

        QuantumMatrix operator()(const double t, const QuantumMatrix &variables) const {
            QuantumMatrix p = variables;
            QuantumMatrix dpdt = -1i * conmutator((*H)(t, p), p)/ hbar;
            for (auto& C : C_ops) {
                QuantumMatrix C_p = (*C)(t, p);
                QuantumMatrix C_p_adj = C_p.adjoint();
                dpdt += 0.5 * (2.0 * C_p * p * C_p_adj - C_p_adj * C_p * p - p * C_p_adj * C_p);
            }
            return dpdt;
        }
    private:
        QobjEvo * H;
        const vector<QobjEvo*> &C_ops;
};


// --------------------------------------------- Solver RK4 -------------------------------------------------
vector<QuantumMatrix> RK4_solver(QobjEvo* derivates, const vector<double> &time_limits, const double step, QuantumMatrix &variable) {
    size_t time_size = (size_t) ((time_limits[1] - time_limits[0]) / step) + 1;
    vector<QuantumMatrix> states;

    states.reserve(time_size);
    states.push_back(variable);
    double t = time_limits[0];

    for (size_t i = 1; i < time_size; i++) {
        QuantumMatrix k1 = step * (*derivates)(t, variable);
        QuantumMatrix k2 = step * (*derivates)(t + 0.5, variable + k1 * 0.5);
        QuantumMatrix k3 = step * (*derivates)(t + 0.5, variable + k2 * 0.5);
        QuantumMatrix k4 = step * (*derivates)(t + step, variable + k3);
        variable += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        states.push_back(variable);
        t += step;
    }
    return states;
}


// --------------------------------------- Master Equation Solver ---------------------------------------------
vector<QuantumMatrix> mesolve(QobjEvo * H, const vector<double> &time_limits, const double step, const QuantumMatrix &p0, const vector<QobjEvo*> &C_ops){
    size_t time_size = (size_t) ((time_limits[1] - time_limits[0]) / step);
    if (time_size < 2){
        cout << "Error: se necesitan al menos dos tiempos para resolver la ecuación" << endl;
        exit(1);
    }
    QuantumMatrix begin_state = p0;

    // Resuelvo la ecuación de Lindblad
    return RK4_solver(new LindbladEvo(H, C_ops), time_limits, step, begin_state);

}


// --------------------------------------- Guardar estados en archivo ---------------------------------------------
void save_states(const vector<QuantumMatrix> &states, const double time_step, const vector<double> &time_limits, const string &filename){
    ofstream file(filename);
    // Guardo primeros los tiempos
    file << "time step: " << time_step << endl;
    file << "time limits: " << time_limits[0] << "," << time_limits[1] << endl;

    //Guardo las matrices de densidades
    file << "states" << endl;
    file << "shape: " << states[0].rows() << "," << states[0].cols() << endl;
    file << "n_states: " << states.size() << endl;

    // 
    for(int i = 0; i < states.size(); i++){
        for(int j = 0; j < states[0].rows(); j++){
            for(int k = 0; k < states[0].cols(); k++){
                if (j == states[0].rows() - 1 && k == states[0].cols() - 1 && i == states.size() - 1)
                    file << states[i](j,k) << endl;
                else
                    file << states[i](j,k) << ",";
            }
        }
    }

    return;

}


// ------------------------------- Creaciones de estados -------------------------------------
QuantumVector basis(const int N, const int i){
    if (N < 0){
        cout << "Error: la dimensión de la base debe ser mayor a 0" << endl;
        exit(1);
    }
    if(i < 0 || i >= N){
        cout << "Error: el índice de la base debe estar entre 0 y N-1" << endl;
        exit(1);
    }
    QuantumVector basis_vector = QuantumVector::Zero(N);

    basis_vector[i] = QuantumValue(1.0, 0.0);

    return basis_vector;
}

QuantumMatrix ket2dm(const QuantumVector &ket){
    if (ket.rows() == 0){
        cout << "Error: el vector de estado no puede ser vacío" << endl;
        exit(1);
    }
    if (ket.cols() == 2){
        cout << "Ya es una DM" << endl;
        return ket; // Si no es un vector columna, se asume que es una matriz de densidad
    }else if (ket.cols() == 1)
        return ket * ket.adjoint();
    else{
        cout << "Error: el vector de estado debe ser un vector columna" << endl;
        exit(1);
    }
}


// Definimos la función de estado de Bell optimizada
QuantumVector bell_state(const std::string &type) {
    QuantumVector bell = QuantumVector::Zero(4);
    const double norm_factor = 1.0 / std::sqrt(2.0);

    if (type == "00") {         // phi_plus
        bell[0] = norm_factor;
        bell[3] = norm_factor;
    } 
    else if (type == "01") {    // phi_minus
        bell[0] = norm_factor;
        bell[3] = -norm_factor;
    } 
    else if (type == "10") {    // psi_plus
        bell[1] = norm_factor;
        bell[2] = norm_factor;
    } 
    else if (type == "11") {    // psi_minus
        bell[1] = norm_factor;
        bell[2] = -norm_factor;
    } 
    else {
        std::cerr << "Error: el estado de Bell no es válido" << std::endl;
        std::exit(1);
    }
    
    return bell;
}
//T. Yu, J.H. Eberly / Optics Communications 264 (2006) 393–397
QuantumMatrix standar_death_form(const vector<double> &diag, const QuantumValue w, const QuantumValue z, const double atol = 1e-6){
    if (diag.size() != 4){
        cout << "Error: la matriz de muerte debe ser de tamaño 4" << endl;
        exit(1);
    }
    double norm = 0;
    for(int i = 0; i < 4; i++)
        norm += diag[i];
    if (abs(norm - 1) > atol){
        cout << "Error: la matriz de muerte debe ser de traza 1: " << norm << endl;
        exit(1);
    }
    QuantumMatrix death_form = QuantumMatrix::Zero(4, 4);

    for(int i = 0; i < 4; i++)
        death_form(i,i) = diag[i];

    death_form(0,3) = w;
    death_form(3,0) = conj(w);

    death_form(1,2) = z;
    death_form(2,1) = conj(z);

    return death_form;
}