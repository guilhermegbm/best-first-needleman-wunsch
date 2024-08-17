#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <pthread.h>

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

// Configurando os valores de forma global para simplificar

//string v = "DRQTAKAAGTD";
//string w = "ERQLAKAAAGTD";
string v = "MVSAIVLYVLLAAAAHSAFASDLDRCTTFDDVQAPNYTQHTSSMRGVYYPDEIFRSDTLYLTQDLFLPFYSNVTGFHTINHTFDNPVIPFKDGIYFAATEKSNVVRGWVFGSTMNNKSQSVIIINNSTNVVIRACNFELCDNPFFAVSKPMGTQTHTMIFDNAFNCTFEYISDAFSLDVSEKSGNFKHLREFVFKNKDGFLYVYKGYQPIDVVRDLPSGFNTLKPIFKLPLGINITNFRAILTAFSPAQDTWGTSAAAYFVGYLKPTTFMLKYDENGTITDAVDCSQNPLAELKCSVKSFEIDKGIYQTSNFRVVPSGDVVRFPNITNLCPFGEVFNATKFPSVYAWERKKISNCVADYSVLYNSTFFSTFKCYGVSATKLNDLCFSNVYADSFVVKGDDVRQIAPGQTGVIADYNYKLPDDFMGCVLAWNTRNIDATSTGNYNYKYRYLRHGKLRPFERDISNVPFSPDGKPCTPPALNCYWPLNDYGFYTTTGIGYQPYRVVVLSFELLNAPATVCGPKLSTDLIKNQCVNFNFNGLTGTGVLTPSSKRFQPFQQFGRDVSDFTDSVRDPKTSEILDISPCSFGGVSVITPGTNASSEVAVLYQDVNCTDVSTAIHADQLTPAWRIYSTGNNVFQTQAGCLIGAEHVDTSYECDIPIGAGICASYHTVSLLRSTSQKSIVAYTMSLGADSSIAYSNNTIAIPTNFSISITTEVMPVSMAKTSVDCNMYICGDSTECANLLLQYGSFCTQLNRALSGIAAEQDRNTREVFAQVKQMYKTPTLKYFGGFNFSQILPDPLKPTKRSFIEDLLFNKVTLADAGFMKQYGECLGDINARDLICAQKFNGLTVLPPLLTDDMIAAYTAALVSGTATAGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKQIANQFNKAISQIQESLTTTSTALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQAAPHGVVFLHVTYVPSQERNFTTAPAICHEGKAYFPREGVFVFNGTSWFITQRNFFSPQIITTDNTFVSGNCDVVIGIINNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKPSGRLVPRGSPGSGYIPEAPRDGQAYVRKDGEWVLLSTFLGHHHHHH";
string w = "MGILPSPGMPALLSLVSLLSVLLMGCVAETGMFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPSRASSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDPPEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELD";

int len_v = v.length();
int len_w = w.length();

vector<vector<double>> matriz_similaridade(len_v + 1, vector<double>(len_w + 1, -numeric_limits<double>::infinity()));
vector<vector<int>> matriz_posicionamento(len_v + 1, vector<int>(len_w + 1, -numeric_limits<double>::infinity()));

MaxHeap heap;

int condicao_parada_alcancada = 0;

// Mutex que garante exclusão mútua para ler e escrever a condição de parada
pthread_mutex_t condicao_parada_mutex;

// Mutex que garante exclusão mútua para leitura e escrita das matrizes de similaridade e posicionamento 
pthread_mutex_t acesso_matrizes_mutex;

// Mutex que garante exclusão mútua ao chamar as funções da heap
pthread_mutex_t heap_mutex;

// Variável de condição que garante, junto à mutex acima, que sempre haverá um nó para a thread processar
pthread_cond_t heap_vazia_cv;

void* parallel_bfnw(void* threadID) {
    long tid = (long) threadID;

    //cout << tid << ": Iniciando" << endl;

    /*
    Não podemos usar a quantidade de elementos na heap (verificar se a heap está vazia) como
    condição de parada, afinal, por exemplo, logo no início da execução, a primeira thread vai
    pegar o primeiro nó para processar e as próximas threads provavelmente vão encontrar uma
    heap vazia, mas isso não significa que o processamento deve ser finalizado - muito pelo
    contrário, acabaram de começar! Portanto, para evitar problemas enquanto não há nós para
    serem processados, devemos criar uma variável de condição que verifica se há elementos na
    heap.
    */
    while (true) {

        // ----- Início seção crítica -----
        // VERIFICANDO a condição de parada, caso outra thread já tenha alcançado o nó sumidouro
        pthread_mutex_lock(&condicao_parada_mutex);
        //cout << tid << ": Entrou SC checando condição de parada" << endl;
        if (condicao_parada_alcancada) {
            // Se alcançou a condição de parada, destrave a mutex ANTES de sair
            //cout << tid << ": Condição de parada ALCANÇADA. Destravando mutex e saindo" << endl;
            pthread_mutex_unlock(&condicao_parada_mutex);
            return NULL;
        }
        //cout << tid << ": Condição de parada NÃO alcançada. Destravando mutex e continuando" << endl;
        pthread_mutex_unlock(&condicao_parada_mutex);
        // Se ainda não alcançou a condição de parada, destrava a mutex e continua
        // ----- Fim seção crítica -----

        // ----- Início seção crítica com Variável de Condição -----
        // Buscando o próximo nó para expandir
        pthread_mutex_lock(&heap_mutex);
        //cout << tid << ": Entrou SC extractMax" << endl;

        while (heap.isEmpty()) { // Boa prática usar "while" em vez de "if"
            //cout << tid << ": Heap vazia! Indo dormir..." << endl;
            pthread_cond_wait(&heap_vazia_cv, &heap_mutex);
            //cout << tid << ": Coisa na Heap (talvez)! Acordando" << endl;
        }

        NodoAuxiliar k = heap.extractMax();

        //cout << tid << ": k encontrado: [" << k.i << ", " << k.j << "]. Destravando mutex heap" << endl;
        pthread_mutex_unlock(&heap_mutex);
        // ----- Fim seção crítica -----

        // ----- Início seção crítica -----
        // Verificando se o valor de k é melhor do que o já calculado
        pthread_mutex_lock(&acesso_matrizes_mutex);
        //cout << tid << ": Entrou SC atualiza matriz" << endl;
        if (k.valor_acumulado > matriz_similaridade[k.i][k.j]) {
            matriz_similaridade[k.i][k.j] = k.valor_acumulado;
            matriz_posicionamento[k.i][k.j] = k.direcao;
            //cout << tid << ": Valor ALTERADO, destrava e segue processando" << endl;
            pthread_mutex_unlock(&acesso_matrizes_mutex);
        } else {
            //cout << tid << ": Nada alterado, destrava e dá continue no loop" << endl;
            pthread_mutex_unlock(&acesso_matrizes_mutex); // Destrava a mutex antes de continuar
            continue;
        }
        // ----- Fim seção crítica -----

        if (k.i == len_v && k.j == len_w) {
            // ----- Início seção crítica -----
            // Sumidouro alcançado! ATUALIZE a condição de parada
            pthread_mutex_lock(&condicao_parada_mutex);
            //cout << tid << ": Entrou SC SUMIDOURO ALCANÇADO!! alterando condição de parada, destravando e saindo" << endl;
            condicao_parada_alcancada = 1;
            pthread_mutex_unlock(&condicao_parada_mutex);
            return NULL;
            // ----- Fim seção crítica -----
        }

        // Processando os novos nós. Essa é a principal parte assíncrona e paralelizável do código
        NodoAuxiliar* no_baixo = nullptr;
        if (k.i < len_v) {
            no_baixo = new NodoAuxiliar(k.i + 1, k.j, k.valor_acumulado - PENALIDADE_INDEL, 0);
        }

        NodoAuxiliar* no_direita = nullptr;
        if (k.j < len_w) {
            no_direita = new NodoAuxiliar(k.i, k.j + 1, k.valor_acumulado - PENALIDADE_INDEL, 1);
        }

        NodoAuxiliar* no_baixo_direita = nullptr;
        if (k.i < len_v && k.j < len_w) {
            int valor_match = matriz_pontuacao_blosum62[dicionario_indice_alfabeto_all_amino.at(v[k.i])][dicionario_indice_alfabeto_all_amino.at(w[k.j])];
            no_baixo_direita = new NodoAuxiliar(k.i + 1, k.j + 1, k.valor_acumulado + valor_match, 2);
        }

        // ----- Início seção crítica -----
        pthread_mutex_lock(&heap_mutex);
        //cout << tid << ": Entrou SC inserção de novos nós na heap" << endl;
        // A inserção na heap exige exclusão mútua para evitar condições de corrida
        if (no_baixo != nullptr) {
            heap.insert(*no_baixo);
        }
        if (no_direita != nullptr) {
            heap.insert(*no_direita);
        }
        if (no_baixo_direita != nullptr) {
            heap.insert(*no_baixo_direita);
        }
        //cout << tid << ": Chamando signal, destravando e recomeçando" << endl;
        // Acordando uma thread que está dormindo na Cond Var da heap vazia
        pthread_cond_signal(&heap_vazia_cv); //TODO: Trocar para broadcast?
        pthread_mutex_unlock(&heap_mutex);
        // ----- Fim seção crítica -----
    }

    //cout << "Fim da thread de id " << tid << endl;
    return NULL;
}

int main(int argc, char* argv[]) {

    int NUM_THREADS;

    if (argc > 1) {
        NUM_THREADS = atoi(argv[1]);
    } else {
        NUM_THREADS = 8;
    }

    condicao_parada_alcancada = 0;
    NodoAuxiliar fonte(0, 0, 0, -numeric_limits<double>::infinity());
    heap.insert(fonte);

    // Inicializando as mutexes, cond vars e atributos
    pthread_mutex_init(&condicao_parada_mutex, NULL);
    pthread_mutex_init(&acesso_matrizes_mutex, NULL);
    pthread_mutex_init(&heap_mutex, NULL);
    pthread_cond_init(&heap_vazia_cv, NULL);

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    pthread_t threads[NUM_THREADS];

    // Criando e inicializando as threads
    for (long tid = 0; tid < NUM_THREADS; tid++) {
        // cout << "Do main: Criando Thread " << tid << endl;
        int return_code = pthread_create(&threads[tid], &attr, *parallel_bfnw, (void*) tid);

        if (return_code) {
            cout << "ERRO; Código de retorno da criação da " << tid << "ª thread foi " << return_code << "!!" << endl;
            return -1;
        }
    }

    // Faz o Join das Threads
    //cout << "Main esperando:" << endl;
    for (long tid = 0; tid < NUM_THREADS; tid++) {
        int joinCode = pthread_join(threads[tid], NULL);
        if (joinCode) {
            cout << "ERRO; Código de retorno da junção da " << tid << "ª thread foi " << joinCode << "!!" << endl;
            return -1;
        }
    }

    cout << "OK" << endl;

    /*cout << "Matriz de Posicionamento:" << endl;
    for (const auto& row : matriz_posicionamento) {
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
    for (const auto& row : matriz_similaridade) {
        for (const auto& elem : row) {
            cout << elem << " ";
        }
        cout << endl;
    }*/

    // Limpando as variáveis
    pthread_attr_destroy(&attr);
    pthread_mutex_destroy(&condicao_parada_mutex);
    pthread_mutex_destroy(&acesso_matrizes_mutex);
    pthread_mutex_destroy(&heap_mutex);
    pthread_cond_destroy(&heap_vazia_cv);

    pthread_exit(NULL);
}
