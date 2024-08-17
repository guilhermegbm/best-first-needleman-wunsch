#include <iostream>
#include <vector>
#include <stdexcept>

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

class MaxHeap {
private:
    vector<NodoAuxiliar> heap;
    int FRONT = 0;

    int parent(int pos) {
        if (pos == 0)
            return 0;
        else
            return (pos - 1) / 2;
    }

    int leftChild(int pos) {
        return (2 * pos) + 1;
    }

    bool hasLeftChild(int pos) {
        return leftChild(pos) < heap.size();
    }

    int rightChild(int pos) {
        return (2 * pos) + 2;
    }

    bool hasRightChild(int pos) {
        return rightChild(pos) < heap.size();
    }

    void heap_swap(int fpos, int spos) {
        swap(heap[fpos], heap[spos]);
    }

    void maxHeapify(int pos) {
        if ((hasLeftChild(pos) && heap[pos] < heap[leftChild(pos)]) || 
            (hasRightChild(pos) && heap[pos] < heap[rightChild(pos)])) {

            if (hasRightChild(pos) && heap[rightChild(pos)] >= heap[leftChild(pos)]) {
                heap_swap(pos, rightChild(pos));
                maxHeapify(rightChild(pos));
            } else {
                heap_swap(pos, leftChild(pos));
                maxHeapify(leftChild(pos));
            }
        }
    }

public:
    int size() {
        return heap.size();
    }

    void insert(NodoAuxiliar element) {
        heap.push_back(element);
        int current = size() - 1;

        while (heap[current] > heap[parent(current)]) {
            heap_swap(current, parent(current));
            current = parent(current);
        }
    }

    bool isEmpty() {
        return heap.size() == 0;
    }

    NodoAuxiliar extractMax() {
        if (isEmpty()) {
            throw runtime_error("Empty Heap");
        }

        NodoAuxiliar first = heap[FRONT];
        NodoAuxiliar last = heap.back();
        heap.pop_back();

        if (!isEmpty()) {
            heap[FRONT] = last;
            maxHeapify(FRONT);
        }

        return first;
    }

    void print() {
        for (const auto& nodo : heap) {
            cout << nodo << endl;
        }
    }
};

int main() {
    MaxHeap maxHeap;

    // Inserindo nodos na heap
    maxHeap.insert(NodoAuxiliar(1, 2, 10, 0));
    maxHeap.insert(NodoAuxiliar(2, 3, 20, 1));
    maxHeap.insert(NodoAuxiliar(3, 4, 5, 2));

    // Imprimindo a heap
    cout << "Max Heap:" << endl;
    maxHeap.print();

    // Extraindo o valor máximo
    cout << "Max Value: " << maxHeap.extractMax() << endl;

    // Imprimindo a heap após a extração
    cout << "Max Heap depois da extração:" << endl;
    maxHeap.print();

    return 0;
}
