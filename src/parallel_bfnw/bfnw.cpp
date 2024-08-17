#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <limits>

using namespace std;

#define NUM_AMINOS 20

#define PENALIDADE_INDEL 0

vector<vector<int>> matriz_pontuacao_blosum62 = {   { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0},
                                                    {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3},
                                                    {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3},
                                                    {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3},
                                                    { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1},
                                                    {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2},
                                                    {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2},
                                                    { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3},
                                                    {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3},
                                                    {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3},
                                                    {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1},
                                                    {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2},
                                                    {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1},
                                                    {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1},
                                                    {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2},
                                                    { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2},
                                                    { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0},
                                                    {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3},
                                                    {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1},
                                                    { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}};

unordered_map<char, int> dicionario_indice_alfabeto_all_amino = {
    {'A', 0},
    {'R', 1},
    {'N', 2},
    {'D', 3},
    {'C', 4},
    {'Q', 5},
    {'E', 6},
    {'G', 7},
    {'H', 8},
    {'I', 9},
    {'L', 10},
    {'K', 11},
    {'M', 12},
    {'F', 13},
    {'P', 14},
    {'S', 15},
    {'T', 16},
    {'W', 17},
    {'Y', 18},
    {'V', 19}
};

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

pair<vector<vector<int>>, vector<vector<double>>> bfnw(
    const string& v, 
    const string& w, 
    const vector<vector<int>>& matriz_pontuacao = matriz_pontuacao_blosum62,
    const unordered_map<char, int>& dicionario_indice_alfabeto = dicionario_indice_alfabeto_all_amino,
    int penalidade_indel = PENALIDADE_INDEL)
{

    int len_v = v.length();
    int len_w = w.length();

    NodoAuxiliar fonte(0, 0, 0, -numeric_limits<double>::infinity());
    MaxHeap heap;
    heap.insert(fonte);

    vector<vector<double>> matriz_similaridade(len_v + 1, vector<double>(len_w + 1, -numeric_limits<double>::infinity()));
    vector<vector<int>> matriz_posicionamento(len_v + 1, vector<int>(len_w + 1, -numeric_limits<double>::infinity()));

    while (!heap.isEmpty()) {
        NodoAuxiliar k = heap.extractMax();

        if (k.valor_acumulado > matriz_similaridade[k.i][k.j]) {
            matriz_similaridade[k.i][k.j] = k.valor_acumulado;
            matriz_posicionamento[k.i][k.j] = k.direcao;
        } else {
            continue;
        }

        if (k.i == len_v && k.j == len_w) {
            break;
        }

        if (k.i < len_v) {
            NodoAuxiliar no_baixo(k.i + 1, k.j, k.valor_acumulado - penalidade_indel, 0);
            heap.insert(no_baixo);
        }

        if (k.j < len_w) {
            NodoAuxiliar no_direita(k.i, k.j + 1, k.valor_acumulado - penalidade_indel, 1);
            heap.insert(no_direita);
        }

        if (k.i < len_v && k.j < len_w) {
            int valor_match = matriz_pontuacao[dicionario_indice_alfabeto.at(v[k.i])][dicionario_indice_alfabeto.at(w[k.j])];
            NodoAuxiliar no_baixo_direita(k.i + 1, k.j + 1, k.valor_acumulado + valor_match, 2);
            heap.insert(no_baixo_direita);
        }
    }

    return {matriz_posicionamento, matriz_similaridade};
}

int main() {
    //string v = "DRQTAKAAGTD";
    //string w = "ERQLAKAAAGTD";
    string v = "MVSAIVLYVLLAAAAHSAFASDLDRCTTFDDVQAPNYTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTINHTFDNPVIPFKDGIYFAATEKSNVVRGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLREFVFKNKDGFLYVYKGYQPIDVVRDLPSGFNTLKPIFKLPLGINITNFRAILTAFSPAQDTWGTSAAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSFEIDKGIYQTSNFRVVPSGDVVRFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSATKLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDSVRDPKTSEILDISPCSFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASYHTVSLLRSTSQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLNRALSGIAAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTTSTALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKPSGRLVPRGSPGSGYIPEAPRDGQAYVRKDGEWVLLSTFLGHHHHHH";
    string w = "MGILPSPGMPALLSLVSLLSVLLMGCVAETGMFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPSRASSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELD";

    auto result = bfnw(v, w);

    cout << "Ok" << endl;
    /*cout << "Matriz de Posicionamento:" << endl;
    for (const auto& row : result.first) {
        for (const auto& elem : row) {
            if (elem == -2147483648) {
                cout << "-inf ";
            } else {
                cout << elem << " ";
            }
            
        }
        cout << endl;
    }

    cout << "Matriz de Similaridade:" << endl;
    for (const auto& row : result.second) {
        for (const auto& elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }*/

    return 0;
}
