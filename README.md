# Lindlbland Solver

Este proyecto contiene un conjunto de herramientas para resolver la ecuación de Lindblad para sistemas cuánticos utilizando el método de Runge-Kutta de cuarto orden (RK4). El código es flexible y permite la simulación de diferentes tipos de sistemas cuánticos, desde un solo qubit hasta sistemas más complejos como qubits entrelazados. El código también incluye funcionalidades para el análisis de procesos como la "muerte súbita del entrelazamiento" (ESD) en sistemas cuánticos acoplados a baños térmicos.

## Descripción

El código resuelve la ecuación de Lindblad numéricamente usando el método de Runge-Kutta de cuarto orden (RK4). Permite simular la evolución de sistemas cuánticos en función de tiempo, considerando efectos de relajación, excitación y desfasaje. El programa puede trabajar con matrices de Pauli y estados de dos niveles (TLS), y permite la creación de matrices de densidad a partir de estados cuánticos representados como vectores columna.

### Características principales:
- Resolución de la ecuación de Lindblad para sistemas cuánticos.
- Implementación del método de Runge-Kutta de cuarto orden (RK4).
- Soporte para matrices de Pauli y estados de Bell.
- Funciones para calcular la evolución temporal y el entrelazamiento de dos qubits.
- Funciones auxiliares para trabajar con productos tensoriales y conmutadores de operadores cuánticos.

## Requisitos

Este proyecto utiliza la biblioteca Eigen para las operaciones matriciales y los cálculos numéricos. Asegúrate de tener instalada la versión más reciente de Eigen.

### Dependencias:
- **Eigen**: Para la manipulación de matrices y vectores. El código utiliza matrices de tipo `MatrixXcd` (matrices de números complejos) y operaciones como el producto tensorial (Kronecker product).
- **C++11 o superior**: Se requiere un compilador que soporte C++11 o versiones posteriores debido al uso de características como `auto` y `nullptr`.

### Instalación de Eigen

Si aún no tienes Eigen, puedes instalarlo de la siguiente manera:

1. Clona el repositorio oficial de Eigen:
   ```bash
   git clone https://gitlab.com/libeigen/eigen.git
   ```
Luego, se puede incluir el directorio eigen en tu proyecto o apuntar a él en tu compilador.

2. Descarga desde la página oficial:
    Alternativamente, puedes descargar Eigen desde su página oficial y seguir las instrucciones de instalación para tu sistema operativo. Después de descargar y descomprimir el archivo, puedes incluir la carpeta Eigen en tu proyecto o apuntar a ella en tu compilador.

### Compilación
Para compilar el código, simplemente usa un compilador C++ que soporte C++11 o superior y asegúrate de que la biblioteca Eigen esté correctamente vinculada.

Ejemplo de compilación usando g++:

    ```bash
    g++ -std=c++11 -I /path/to/eigen your_file.cpp -o quantum_solver
    ```

### Funciones principales

 - ```RK4_solver```: Resuelve la ecuación de Lindblad utilizando el método de Runge-Kutta de cuarto orden.
 - ```mesolve```: Resolver la ecuación de Lindblad para un sistema cuántico general, dado un operador Hamiltoniano y operadores de colapso.
- ```save_states```: Guarda los estados evolutivos del sistema en un archivo de texto.
- ```basis```: Genera el vector base correspondiente a un índice específico de un espacio de Hilbert de dimensión N.
- ```ket2dm```: Convierte un estado cuántico en un vector columna (ket) a una matriz de densidad (dm).
- ```bell_state```: Genera estados de Bell (entrelazados) con diferentes tipos.
- ```standar_death_form```: Genera una matriz de proceso para estudiar la muerte súbita del entrelazamiento según Yu y Eberly (2006).

#### Autor

Este proyecto fue creado por Maximiliano Gatto, como parte de un estudio en la simulación numérica de sistemas cuánticos.