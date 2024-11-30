#include <iostream>
#include <vector>
// #include <utils.cpp>
#include <complex>
#include <fstream>

using namespace std;


// Definiciones más fisicas
using quantum_value = complex<double>;
using quantum_object = vector<vector<quantum_value>>;

// Defino matrices de Pauli
quantum_object sigma_x = {{0, 1}, {1, 0}};
quantum_object sigma_y = {{0, -1i}, {1i, 0}};
quantum_object sigma_z = {{1, 0}, {0, -1}};
quantum_object sigma_p = {{0, 1}, {0, 0}};
quantum_object sigma_m = {{0, 0}, {1, 0}};
quantum_object sigma_0 = {{1, 0}, {0, 1}};
quantum_object I = {{1, 0}, {0, 1}};

// Multiplicación de matrices
quantum_object matrix_multiply(const quantum_object &q1, const quantum_object &q2){
    vector<size_t> shape1 = {q1.size(), q1[0].size()};
    vector<size_t> shape2 = {q2.size(), q2[0].size()};

    if (shape1[1] != shape2[0]){
        cout << "Error: las matrices no son compatibles para multiplicar" << endl;
        exit(1);
    }

    quantum_object result(shape1[0], vector<quantum_value>(shape2[1], 0));
    for(int i = 0; i < shape1[0]; i++)
        for(int j = 0; j < shape2[1]; j++)
            for(int k = 0; k < shape1[1]; k++)
                result[i][j] += q1[i][k] * q2[k][j];
    
    return result;
}

// Definición de estados
quantum_object basis(const int N, const int i){
    if (N < 0){
        cout << "Error: la dimensión de la base debe ser mayor a 0" << endl;
        exit(1);
    }
    if(i < 0 || i >= N){
        cout << "Error: el índice de la base debe estar entre 0 y N-1" << endl;
        exit(1);
    }
    quantum_object basis(N, vector<quantum_value>(1, 0));
    basis[i][0] = 1;
    return basis;
}

// daggar un operador o estado (cambia entre bra y ket)
quantum_object dag(const quantum_object &q){
    vector<size_t> shape = {q.size(), q[0].size()};
    quantum_object result(shape[1], vector<quantum_value>(shape[0], 0));

    for(int i = 0; i < shape[0]; i++)
        for(int j = 0; j < shape[1]; j++)
            result[j][i] = conj(q[i][j]);
    
    return result;
}


vector<vector<double>> abs(const quantum_object &q){
    vector<size_t> shape = {q.size(), q[0].size()};
    vector<vector<double>> result(shape[0], vector<double>(shape[1], 0));

    for(int i = 0; i < shape[0]; i++)
        for(int j = 0; j < shape[1]; j++)
            result[i][j] = abs(q[i][j]);
    
    return result;
}

// Defino operaciones
quantum_object operator*(const quantum_object &q1, const quantum_object &q2){
    return matrix_multiply(q1, q2);
}

quantum_object conmutator(const quantum_object &q1, const quantum_object &q2){
    return q1 * q2 - q2 * q1;
}

quantum_object operator+=(quantum_object &q1, const quantum_object &q2){
    vector<size_t> shape1 = {q1.size(), q1[0].size()};
    vector<size_t> shape2 = {q2.size(), q2[0].size()};

    if (shape1[0] != shape2[0] || shape1[1] != shape2[1]){
        cout << "Error: las matrices deben ser del mismo tamaño para sumar" << endl;
        exit(1);
    }

    for(int i = 0; i < shape1[0]; i++)
        for(int j = 0; j < shape1[1]; j++)
            q1[i][j] += q2[i][j];
    
    return q1;
}

quantum_object ket2dm(const quantum_object &q){
    if (q[0].size() != 1)   // Si no es un ket, no se puede convertir a densidad
        return q;
        
    quantum_object result = q * dag(q);
    return result;
}

quantum_object normalize(const quantum_object &q1){
    quantum_object result = q1;
    quantum_value norm = 0;
    for(int i = 0; i < q1.size(); i++)
        for(int j = 0; j < q1[0].size(); j++)
            norm += q1[i][j] * conj(q1[i][j]);
    
    norm = sqrt(norm);
    for(int i = 0; i < q1.size(); i++)
        for(int j = 0; j < q1[0].size(); j++)
            result[i][j] = q1[i][j]/norm;
    
    return result;
}

bool is_normalized(const quantum_object &q, double atol = 1e-10){
    double norm = 0;
    cout << "normalizing.." << endl;
    for(int i = 0; i < q.size(); i++)
        for(int j = 0; j < q[0].size(); j++)
            norm += abs(q[i][j]);

    norm = sqrt(norm);
    cout << "norm: " << norm << endl;
    if (abs(norm - 1) < atol)
        return true;
    else
        return false;
}

quantum_object tensor(const quantum_object &q1, const quantum_object &q2){
    vector<size_t> shape1 = {q1.size(), q1[0].size()};
    vector<size_t> shape2 = {q2.size(), q2[0].size()};

    quantum_object result(shape1[0]*shape2[0], vector<quantum_value>(shape1[1]*shape2[1], 0));

    for(int i = 0; i < shape1[0]; i++)
        for(int j = 0; j < shape1[1]; j++)
            for(int k = 0; k < shape2[0]; k++)
                for(int l = 0; l < shape2[1]; l++)
                    result[i*shape2[0] + k][j*shape2[1] + l] = q1[i][j] * q2[k][l];
    
    return result;
}


//------------  Clase base abstracta para los functors ---------------
class QobjEvo {
    public:
        virtual quantum_object operator()(const double t, const quantum_object &variables) const = 0;
};

class LindbladEvo : public QobjEvo {
    public:
        LindbladEvo(QobjEvo * H, const vector<QobjEvo*> &C_ops) : H(H), C_ops(C_ops) {}

        quantum_object operator()(const double t, const quantum_object &variables) const {
            quantum_object p = variables;
            quantum_object dpdt = -1i * conmutator((*H)(t, p), p)/ hbar;
            for(int i = 0; i < C_ops.size(); i++){
                dpdt += 0.5 * (2.0 * (*C_ops[i])(t, p) * p * dag((*C_ops[i])(t, p)) - dag((*C_ops[i])(t, p)) * (*C_ops[i])(t, p) * p - p * dag((*C_ops[i])(t, p)) * (*C_ops[i])(t, p));
            }
            return dpdt;
        }
    private:
        QobjEvo * H;
        vector<QobjEvo*> C_ops;
};

// ----------------------------- Solver de RK4 ----------------------------
quantum_object RK4_solver(QobjEvo* derivates, double t, quantum_object &variable, double h) {
    // verifico que las derivadas sean del mismo tamaño que las variables
    const vector<size_t> shape = {variable.size(), variable[0].size()}; // shape de la density matrix

    quantum_object k1(shape[0], vector<quantum_value>(shape[1], 0)), k2(shape[0], vector<quantum_value>(shape[1], 0)), k3(shape[0], vector<quantum_value>(shape[1], 0)), k4(shape[0], vector<quantum_value>(shape[1], 0));

    // Calculo los k1
    k1 = h * (*derivates)(t, variable);
    // Calculo los k2
    k2 = h * (*derivates)(t + h*0.5, variable + k1*0.5);
    // Calculo los k3
    k3 = h * (*derivates)(t + h*0.5, variable + k2*0.5);
    // Calculo los k4
    k4 = h * (*derivates)(t + h, variable + k3);

    return variable + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;

}


//Master Equation Solve
vector<quantum_object> mesolve(QobjEvo * H, const vector<double> &times, const quantum_object &p0, const vector<QobjEvo*> &C_ops){
    size_t time_size = times.size();
    if (time_size < 2){
        cout << "Error: se necesitan al menos dos tiempos para resolver la ecuación" << endl;
        exit(1);
    }
    quantum_object begin_state = ket2dm(p0);

    vector<quantum_object> states(time_size, begin_state);   // Estados a resolver

    // Defino las derivadas como un QobjEvo
    LindbladEvo * lindbladian = new LindbladEvo(H, C_ops); 

    // Resuelvo la ecuación de Lindblad
    for(int i = 1; i < time_size; i++){
        states[i+1] = RK4_solver(lindbladian, times[i], states[i], times[i+1] - times[i]);
    }

    delete lindbladian;

    return states;
}

void save_states(const vector<quantum_object> &states, const double time_step, const vector<double> &time_limits, const string &filename){
    ofstream file(filename);
    // Guardo primeros los tiempos
    file << "time step: " << time_step << endl;
    file << "time limits: " << time_limits[0] << "," << time_limits[1] << endl;

    //Guardo las matrices de densidades
    file << "states" << endl;
    file << "shape: " << states[0].size() << "," << states[0][0].size() << endl;
    file << "n_states: " << states.size() << endl;

    // 
    for(int i = 0; i < states.size(); i++){
        for(int j = 0; j < states[0].size(); j++){
            for(int k = 0; k < states[0][0].size(); k++){
                if (j == states[0].size() - 1 && k == states[0][0].size() - 1 && i == states.size() - 1)
                    file << states[i][j][k] << endl;
                else
                    file << states[i][j][k] << ",";
            }
        }
    }

    return;

}
