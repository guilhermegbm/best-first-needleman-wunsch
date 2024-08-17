import numpy as np

import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

import utils.constantes as constantes

def needleman_wunsch_iterativo(
    v,
    w,
    matriz_pontuacao=constantes.matriz_pontuacao_blosum62,
    dicionario_indice_alfabeto=constantes.dicionario_indice_alfabeto_all_amino, 
    penalidade_indel=constantes.PENALIDADE_INDEL):
    """
    Implementação iterativa do Algoritmo Needleman-Wunsch para construir alinhamentos
    globais par-a-par entre sequências.
    """

    len_v = len(v)
    len_w = len(w)

    # A matriz que armazena as regras de posicionamento ou direcionamento do caminho
    # 0 = "north", 1 = "west" e 2 = "northwest"
    matriz_posicionamento = np.zeros((len_v+1, len_w+1))
    matriz_posicionamento[0, :] = np.ones((len_w+1))


    # Inicializando a matriz de recorrencia ou memoização
    matriz_similaridade = np.zeros((len_v+1, len_w+1))
    for i in range(1, len_v+1):
        matriz_similaridade[i, 0] = matriz_similaridade[i-1, 0] - penalidade_indel

    for j in range(1, len_w+1):
        matriz_similaridade[0, j] = matriz_similaridade[0, j-1] - penalidade_indel

    for i in range(1, len_v+1):
        for j in range(1, len_w+1):

            # Inserção em v:
            valor_insercao_v = matriz_similaridade[i-1, j] - penalidade_indel

            # Inserção em w:
            valor_insercao_w = matriz_similaridade[i, j-1] - penalidade_indel

            # Casamento (Se vai dar match mesmo ('T' == 'T') ou deu mismatch ('G' != 'T'),
            # não importa. Isso está mapeado na matriz de pontuação)
            valor_match = matriz_similaridade[i-1, j-1] + matriz_pontuacao[dicionario_indice_alfabeto[v[i-1]], dicionario_indice_alfabeto[w[j-1]]]

            # Atualizando o valor na matriz de similaridade
            matriz_similaridade[i, j] = max(valor_insercao_v, valor_insercao_w, valor_match)

            # Atualizando a operação que deve ser feita no caminhamento
            if matriz_similaridade[i, j] == valor_match:
                matriz_posicionamento[i, j] = 2
            elif matriz_similaridade[i, j] == valor_insercao_w:
                matriz_posicionamento[i, j] = 1
            elif matriz_similaridade[i, j] == valor_insercao_v:
                matriz_posicionamento[i, j] = 0

            """
            Observação IMPORTANTE: A ordem dos if-elif acima INFLUENCIAM no caminhamento final, no sentido de
            que alguns caminhos terão prioridade quando houver algum empate nos valores. Por exemplo, para o
            exemplo no livro (Jones, Pevzner, p. 173), a configuração acima dá um resultado DIFERENTE (mas eu
            acho que é tão bom quanto!), porém, ao trocar a ordem de verificação do valor_insercao_w pelo
            valor_insercao_v acima, o resultado se iguala ao exemplo do livro!
            """

    return matriz_posicionamento, matriz_similaridade