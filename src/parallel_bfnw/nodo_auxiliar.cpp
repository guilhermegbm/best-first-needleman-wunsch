#include <iostream>
#include <string>

using namespace std;

struct NodoAuxiliar {
    int i;
    int j;
    int valor_acumulado;
    int direcao;

    // Construtor
    NodoAuxiliar(int i_, int j_, int valor_acumulado_, int direcao_)
        : i(i_), j(j_), valor_acumulado(valor_acumulado_), direcao(direcao_) {}

    // Sobrecarga do operador de inserção para imprimir a struct
    friend ostream& operator<<(ostream& os, const NodoAuxiliar& nodo) {
        os << "i: " << nodo.i << ", j: " << nodo.j << ", va: " << nodo.valor_acumulado << ", dir: " << nodo.direcao;
        return os;
    }

    // Sobrecarga dos operadores de comparação
    bool operator<(const NodoAuxiliar& other) const {
        return valor_acumulado < other.valor_acumulado;
    }

    bool operator<=(const NodoAuxiliar& other) const {
        return valor_acumulado <= other.valor_acumulado;
    }

    bool operator>(const NodoAuxiliar& other) const {
        return valor_acumulado > other.valor_acumulado;
    }

    bool operator>=(const NodoAuxiliar& other) const {
        return valor_acumulado >= other.valor_acumulado;
    }
};

int main() {
    // Exemplo de uso
    NodoAuxiliar nodo1(1, 2, 10, 0);
    NodoAuxiliar nodo2(2, 3, 20, 1);

    cout << nodo1 << endl;
    cout << nodo2 << endl;

    if (nodo1 < nodo2) {
        cout << "Nodo1 tem valor acumulado menor que Nodo2" << endl;
    }

    if (nodo1 <= nodo2) {
        cout << "Nodo1 tem valor acumulado menor ou igual a Nodo2" << endl;
    }

    if (nodo1 > nodo2) {
        cout << "Nodo1 tem valor acumulado maior que Nodo2" << endl;
    }

    if (nodo1 >= nodo2) {
        cout << "Nodo1 tem valor acumulado maior ou igual a Nodo2" << endl;
    }

    return 0;
}
