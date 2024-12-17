import copy

import numpy as np


def transpose(A):
    # transposta da matriz A
    A_t = []
    A_t_linha = []
    if type(A) == "int":
        m_A = 1
    else:
        m_A = len(A)
    if type(A[0]) == "int":
        n_A = 1
    else:
        n_A = len(A[0])

    for j in range(n_A):
        A_t_linha = []
        for i in range(m_A):
            A_t_linha.append(A[i][j])
        A_t.append(A_t_linha)
    return A_t


def vectorProduct(A, B):
    # produto do vetor linha A e do vetor coluna B, o que resulta em um número
    # Note, que A vai ser uma lista normal, enquanto B vai ser uma lista de listas
    m_A = 1
    n_A = len(A)
    m_B = len(B)
    n_B = 1
    c = 0
    if n_A == m_B:
        for i in range(n_A):
            c += A[i] * B[i][0]
    else:
        c = "Tamanhos incompatíveis"
    return c


def verificafase1(A, b):
    m = len(A)
    n = len(A[0])
    q = n - m
    r = n - q
    count_lines = 0
    numbers_line = 0
    indexes = []
    for i in range(m):
        # reseta o contador de números positivos na linha
        numbers_line = 0
        for j in range(r):
            # verifica quantos não zeros tem na linha e em quais posições (algo semelhante a uma identidade)
            if A[i][j+q] > 0:
                if j+q not in indexes:
                    numbers_line += 1
                    indexes.append(j+q)
        # primeira linha só tem um número nas n-m colunas finais
        if numbers_line == 1:
            count_lines += 1
    # se todas as linhas tem exatamente um elemento não-nulo
    if count_lines == r:
        # gerando a partição básica e listas para armazenar seus índices
        B = []
        N = []
        indB = []
        indN = []
        for i in range(m):
            B.append(A[i][q:])
            N.append(A[i][:q])
        for j in range(r):
            indB.append(j+q)
        for k in range(q):
            indN.append(k)
        # resolvendo o sistema linear para encontrar a solução básica
        x = np.linalg.solve(B, b)
        # verificando se a solução básica é factível
        for j in range(len(x)):
            if x[j] < 0:
                return 0, 0, 0, 0
        return B, N, indB, indN
    else:
        return 0, 0, 0, 0


def solucaobasica(B, b):
    # encontra a solução básica da partição básica atual
    x = np.linalg.solve(B, b)
    for j in range(len(x)):
        if x[j] < 0:
            return 0
    return x


def vetormultiplicador(B, c_B):
    B_t = transpose(B)
    Lambda = np.linalg.solve(B_t, c_B)
    return Lambda


def custosRelativos(c_N, Lambda, N):
    # transpondo N, para ficar mais fácil de trabalhar
    N_t = transpose(N)
    # lista com os custos relativos
    c_r = []
    for i in range(len(c_N)):
        c_r_elem = c_N[i][0] - vectorProduct(N_t[i], Lambda)
        c_r.append(c_r_elem)
    return c_r


def direcaoSimplex(B, a_N):
    d = np.linalg.solve(B, a_N)
    return d


def tamanhoPasso(d, x_B):
    # encontrando o tamanho do passo e quem sai da base. epsilon_ind guarda as posições dos elementos do vetor direção positivos
    epsilon_ind = []
    # epsilon_list guarda os epsilons possíveis
    epsilon_list = []
    for i in range(len(d)):
        if d[i] > 0:
            epsilon_ind.append(i)
    if epsilon_ind == []:
        # nesse caso, nenhum elemento do vetor diretor é positivo, então a função retorna 2 números negativos
        # para indicar que o problema é ilimitado
        return -42, -42
    else:
        for j in epsilon_ind:
            epsilon_list.append(x_B[j]/d[j])
        epsilon = min(epsilon_list)
        x_sai_ind = epsilon_ind[epsilon_list.index(epsilon)]
        return epsilon, x_sai_ind


def simplex(A, b, c, B, N, indB, indN):
    # variável contadora de interações
    inte = 1
    while True:
        # os parâmetros do problema e uma solução básica factível inicial, além das colunas da partição
        c_B = []
        c_N = []
        for j in indB:
            c_B.append(c[j])
        for k in indN:
            c_N.append(c[k])
        x_B = solucaobasica(B, b)
        Lambda = vetormultiplicador(B, c_B)
        c_r = custosRelativos(c_N, Lambda, N)
        # custo relativo mínimo
        c_N_k = min(c_r)
        if c_N_k >= 0:
            return x_B, indB, indN
        else:
            # encontrando o índice que entrará na base na próxima iteração
            var_N = 0
            for j in range(len(c_r)):
                if c_N_k == c_r[j]:
                    var_N = indN[j]
        # pegando a coluna a_N_k
        A_t = transpose(A)
        a_N_k = A_t[var_N]
        # vetor diretor simplex: d (y será para variáveis artificiais)
        d = direcaoSimplex(B, a_N_k)
        # tamanho do base e índice da variável básica que sairá da base
        epsilon, x_sai_ind = tamanhoPasso(d, x_B)
        if epsilon < 0:
            return 0, 0, 0
        # atualização da base
        else:
            # atualizando as listas com os índices básicos e não-básicos
            indB_entra = var_N
            indB_sai = indB[x_sai_ind]
            indB.remove(indB_sai)
            indB.append(indB_entra)
            indN.remove(indB_entra)
            indN.append(indB_sai)
            B = []
            N = []
            B_t = []
            N_t = []
            for i in range(len(indB)):
                B_t.append(A_t[indB[i]])
            B = transpose(B_t)
            for j in range(len(indN)):
                N_t.append(A_t[indN[j]])
            N = transpose(N_t)
            inte += 1


def funcaoObjetivo(X, c):
    f = vectorProduct(X, c)
    return f


#Input para a matriz A#
Rest = int(input("Número de restrições: "))

Vari = int(input("Número de variáveis: "))
print("--------------------------------------")
print("informe os coeficientes da matriz A:")
A = []


for i in range(Rest):

    MatJunt = list(map(float, input().split()))

    A.append(MatJunt)
#Input para o vetor b#

b = []
Restb = Rest
Varib = 1
print("--------------------------------------")
print("informe os coeficientes do vetor b: ")

for i in range(Restb):
    MatJuntb = list(map(float, input().split()))
    b.append(MatJuntb)


for j in range(len(b)):
    if b[j][0] < 0:
        b[j][0] = -b[j][0]
        for i in range(len(A[0])):
            A[j][i] = -A[j][i]
#Input para o custo#
c = []
Restc = Vari
Varic = Vari
print("--------------------------------------")
print("informe os coeficientes da Função Objetivo: ")

for i in range(Restc):
    MatJuntc = list(map(float, input().split()))
    c.append(MatJuntc)
print("--------------------------------------")
m = len(A)  # número de linhas
n = len(A[0])  # número de colunas
B, N, indB, indN = verificafase1(A, b)
if B != 0:
    x_B, indB, indN = simplex(A, b, c, B, N, indB, indN)
    if indB == 0:
        print("O problema em questão é ilimitado")
    else:
        # criando um vetor X para armazenar a resposta
        X = []
        for i in range(n):
            X.append(0)
        for j in range(len(indB)):
            X[indB[j]] = x_B[j][0]
        f = funcaoObjetivo(X, c)
        print("-------------------------------")
        print("Solução Básica factível ótima:")
        print(X)
        print("-------------------------------")
        print("Valor da solução:")
        print(f)
        print("-------------------------------")
else:
    # método das duas fases
    c_aux = []
    A_aux = copy.deepcopy(A)
    for i in range(len(c)):
        c_aux.append([0])
    # criando as variáveis artificiais
    for j in range(m):
        for i in range(m):
            A_aux[j].append(0)
        c_aux.append([1])
    for i in range(m):
        A_aux[i][i+n] = 1
    print(A_aux)
    A_aux_t = transpose(A_aux)
    # partição básica inicial para o método das duas fases
    B_duas, N_duas, indB_duas, indN_duas = verificafase1(A_aux, b)
    # executa o método simplex para o problema auxiliar
    x_B_duas, indB_duas, indN_duas = simplex(
        A_aux, b, c_aux, B_duas, N_duas, indB_duas, indN_duas)
    # criando um vetor X_duas para armazenar a resposta do método das duas fases
    X_duas = []
    for i in range(len(A_aux[0])):
        X_duas.append(0)
    for j in range(len(indB_duas)):
        X_duas[indB_duas[j]] = x_B_duas[j][0]
    g = funcaoObjetivo(X_duas, c_aux)
    # verifica se o valor ótimo é diferente de 0
    if g == 0:
        lista_aux = []
        for j in range(len(indN_duas)):
            if indN_duas[j] >= len(A[0]):
                lista_aux.append(indN_duas[j])
        for elem in lista_aux:
            indN_duas.remove(elem)
        # criando a partição básica factível encontrada pelo método das duas fases
        B = []
        N = []
        B_t = []
        N_t = []
        A_t = transpose(A)
        for i in range(len(indB_duas)):
            B_t.append(A_t[indB_duas[i]])
        B = transpose(B_t)
        for j in range(len(indN_duas)):
            N_t.append(A_t[indN_duas[j]])
        N = transpose(N_t)
        # executando o método simplex para o problema original a partir da partição básica encontrada pelo método
        x_B, indB, indN = simplex(A, b, c, B, N, indB_duas, indN_duas)
        # criando um vetor X para armazenar a resposta
        if indB == 0:
            print("O problema em questão é ilimitado")
        else:
            # criando um vetor X para armazenar a resposta
            X = []
            for i in range(n):
                X.append(0)
            for j in range(len(indB)):
                X[indB[j]] = x_B[j][0]
            f = funcaoObjetivo(X, c)
            print("-------------------------------")
            print("Solução Básica factível ótima:")
            print(X)
            print("-------------------------------")
            print("Valor da solução:")
            print(f)
            print("-------------------------------")
    else:
        print("Problema infactível")