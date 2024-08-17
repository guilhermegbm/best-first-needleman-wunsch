#include <pthread.h>
#include <iostream>
#include <string>

using namespace std;

#define NUM_THREADS 5

void* printar_hello(void* threadID) {
    long tid = (long) threadID;

    cout << "Hello World da Thread " << tid << endl;
    pthread_exit(NULL);
}

int main (int argc, char *argv[]) {
    pthread_t threads[NUM_THREADS];
    int rc;
    
    for (long t = 0; t < NUM_THREADS; t++) {
        cout << "Do main: Criando Thread " << t << endl;
        rc = pthread_create(&threads[t], NULL, *printar_hello, (void *)t);

        if (rc) {
            cout << "ERRO; Código de retorno da criação da "<< t << "ª thread foi " << rc << " !!" << endl;
            //exit(-1);
            return -1;
        }
    }

    pthread_exit(NULL);
}