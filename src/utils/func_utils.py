def computar_mapa_posicionamento(matriz_posicionamento):
    """
    Dada a matriz bidimensional de posicionamento gerada após um algoritmo de
    um alinhamento par-a-par, essa função printa no console o mapa de posicionamento
    que descreve os passos para realizar o caminhamento para cada elemento 
    processado na matriz. Os valores possíveis nessa matriz são 2 (direção
    "northwest"), 1 ("west") e 0 ("north"). Para qualquer outro valor (como
    np.nan, -np.inf, etc), assume-se que o elemento não foi processado e um
    espaço em branco " " será printado.
    """
    
    str_mapa = ""
    for i in range(matriz_posicionamento.shape[0]):
        for j in range(matriz_posicionamento.shape[1]):
            if matriz_posicionamento[i, j] == 2:
                str_mapa += "\\ "
            elif matriz_posicionamento[i, j] == 1:
                str_mapa += "_ "
            elif matriz_posicionamento[i, j] == 0:
                str_mapa += "| "
            else:
                str_mapa += "  " # Isso inclui -inf

        str_mapa+="\n"

    return str_mapa

def computar_alinhamento(v, w, matriz_posicionamento):
    """
    Dadas as sequências v e w originais e a matriz de posicionamento gerada após
    o alinhamento de v e w, essa função computa a representação do alinhamento
    par-a-par entre as sequências, adicionando "_" em gaps.
    """
    str_v = ""
    str_w = ""

    contador_v = len(v)
    contador_w = len(w)

    while contador_v > 0 or contador_w > 0:
        if matriz_posicionamento[contador_v, contador_w] == 2:
            str_v += v[contador_v-1]
            str_w += w[contador_w-1]
            contador_v-=1
            contador_w-=1
        elif matriz_posicionamento[contador_v, contador_w] == 1:
            str_v += "-"
            str_w += w[contador_w-1]
            contador_w-=1
        elif matriz_posicionamento[contador_v, contador_w] == 0:
            str_v += v[contador_v-1]
            str_w += "-"
            contador_v-=1

    return str_v[::-1], str_w[::-1]

import random

"""def gerar_string_aleatoria(tam, dicionario):
    chaves = list(dicionario.keys())
    string = ""

    for i in range(tam):
        string += random.choice(chaves)

    return string"""

def gerar_string_aleatoria(tam, dicionario, seed=None):
    """
    Dado um tamanho "tam" e um dicionário de caracteres, essa função gera uma
    string aleatória de tamanho tam com os caracteres que podem ser encontrados
    no dicionário. Opcionalmente, passe um parâmetro numérico "seed" para
    configurar a semente do gerador aleatório
    """
    if seed != None:
        random.seed(seed)
    
    chaves = list(dicionario.keys())
    string = [""] * tam

    for i in range(tam):
        string[i] = random.choice(chaves)

    return "".join(string)

def mutar_string(string_original, dicionario, k):
    """
    Dada uma string (chamada aqui de "string_original"), um dicionário contendo
    caracteres de um alfabeto e um valor numérico 0 <= k <= 100, essa função
    aplica uma mutação na string_original, gerando uma nova string com k% de
    seus elementos alterados para caracteres amostrados aleatoriamente do
    dicionário. As k% posições mutadas também são aleatórias.
    """

    qtde_posicoes_para_mutar = int(len(string_original) * (k / 100))

    len_string_original = len(string_original)

    string_mutada = list(string_original) # A string mutada começa idêntica à string original

    chaves = list(dicionario.keys())

    for k in range(0, qtde_posicoes_para_mutar):
        index_aleatorio = random.randint(0, len_string_original-1)
        # NOTA 1: A mesma posição aleatória pode ser gerada (e mutada) mais de uma vez!

        string_mutada[index_aleatorio] = random.choice(chaves)
        # NOTA 2: Um valor pode ser aleatoriamente mutado para ele mesmo!

    return "".join(string_mutada)

def computar_porcentagem_identidade(str_orig, str_mut):
    """
    Dadas duas sequências alinhadas (ou, no mínimo, de mesmo tamanho), essa função computa
    a métrica de identidade entre as strings, que corresponde à porcentagem de similaridade
    entre seus elementos, ou seja, quantos por cento dos seus caracteres na mesma posição
    são iguais
    """
    str_len = len(str_orig)

    num_ident = 0
    for i in range(str_len):
        if str_orig[i] == str_mut[i]:
            num_ident+=1

    return (num_ident / str_len) * 100