#include <vector>
#include <fstream>
#include <string>
#define hbar 1.0

using namespace std;

// ---------------- Sobrecarga de operadores --------------------
template <typename T>
vector<T> operator*(const vector<T> &a, const vector<T> &b) {
    const size_t size = a.size();
if (size != b.size()){
        string error_msg = "Error: No se puede multiplicar un vector de longitud " + to_string(size) + " con un vector de longitud " + to_string(b.size()) + ".";

        throw std::invalid_argument(error_msg);
    }
    vector<T> result(size);  // Define el tamaño directamente
    for (size_t i = 0; i < size; ++i)
        result[i] = a[i] * b[i];

    return result;
}

template <typename T>
vector<T> operator/(const vector<T> &a, const vector<T> &b) {
    const size_t size = a.size();
    if (size != b.size()){
        string error_msg = "Error: No se puede dividir un vector de longitud " + to_string(size) + " con un vector de longitud " + to_string(b.size()) + ".";
        throw std::invalid_argument(error_msg);
    }
    
    vector<T> result(size);
    for(int i = 0; i < size; i++)
        result[i] = a[i] / b[i];
    
    return result;
}

template <typename T>
vector<T> operator+(const vector<T> &a, const vector<T> &b) {
    const size_t size = a.size();
    if (size != b.size()){
        string error_msg = "Error: No se puede sumar un vector de longitud " + to_string(size) + " con un vector de longitud " + to_string(b.size()) + ".";
        throw std::invalid_argument(error_msg);
    }
    vector<T> result(size);
    for(int i = 0; i < size; i++)
        result[i] = a[i] + b[i];

    return result;
}

template <typename T>
vector<T> operator-(const vector<T> &a, const vector<T> &b) {
    const size_t size = a.size();
    if (size != b.size()){
        string error_msg = "Error: No se puede restar un vector de longitud " + to_string(size) + " con un vector de longitud " + to_string(b.size()) + ".";

        throw std::invalid_argument(error_msg);
    }
    vector<T> result(size);
    for(int i = 0; i < a.size(); i++) 
        result[i] = a[i] - b[i];

    return result;
}

template <typename T, typename U>
vector<T>  operator*(const U &c, const vector<T> &a) {
    const size_t size = a.size();
    vector<T> result(size);
    for(int i = 0; i < a.size(); i++)
        result[i] = c * a[i];

    return result;
}

template <typename T, typename U>
vector<T> operator*(const vector<T> &a, const U &c) {
    const size_t size = a.size();
    vector<T> result(size);
    for(int i = 0; i < a.size(); i++)
        result[i] = c * a[i];

    return result;
}

template <typename T, typename U>
vector<T> operator/(const vector<T> &a, const U &c) {
    const size_t size = a.size();
    vector<T> result(size);
    for(int i = 0; i < a.size(); i++)
        result[i] = a[i]/c;

    return result;
}

template <typename T, typename U>
vector<T> operator/(const U &c, const vector<T> &a) {
    const size_t size = a.size();
    vector<T> result(size);
    for(int i = 0; i < a.size(); i++)
        result[i] = a[i]/c;

    return result;
}

template <typename T>
basic_ostream<char> & operator<<(basic_ostream<char> &os, const vector<T> &v) {
    for(size_t i = 0; i < v.size(); i++) {
        if (i == v.size() - 1) 
            os << v[i] << endl;
        else
            os << v[i] << ',';
    }
    return os;
}
