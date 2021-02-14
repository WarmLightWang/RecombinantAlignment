import sys
import submatrix


def score_recombine_align(alignment, parents, substitution_matrix, r):
    """Scores a recombination alignment.

    Args:
        alignment: An alignment of three sequences: x, y_0, and y_1 represented as a list of three strings
        parents: a string of '0' and '1' indicating the parent
                 sequence (y_0 or y_1) at each column in the alignment.
        substitution_matrix: the substitution matrix, represented as a dictionary 
        r: the recombination penalty (generally <= 0)
    Returns:
        The score (a number) of alignment given the substitution matrix and r
    """
    score = 0
    x_aligned = alignment[0]
    y_aligned = alignment[1:3]
    parent_index = list(map(int, parents))
    for i in range(len(x_aligned)):
        score += substitution_matrix[(x_aligned[i],
                                      y_aligned[parent_index[i]][i])]
        if i > 0 and parent_index[i] != parent_index[i - 1]:
            score += r
    return score


def recombine_align(x, y, substitution_matrix, r):
    """Computes an optimal global pairwise recombination alignment of sequence x
       with a pair of aligned parental sequences y

    Uses a linear gap scoring function with the space score encoded in the given substitution matrix.
    In the case of multiple optimal alignments, the order of preference for the traceback is given
    by the order of the cases in the dynamic programming recurrences.

    Args:
        x: a string
        y: alignment of two strings, representing the possible parental sequences of x.
        substitution_matrix: the substitution matrix, represented as a dictionary
        r: the recombination penalty (generally <= 0)
    Returns:
        A tuple of the form (score, alignment, parents), where score is the optimal score,
        alignment is a multiple alignment of x and y (represented as a list of three aligned sequences),
        and parents is a string consisting of the characters '0' and '1' with length equal to the number
        of columns of the multiple alignment, indicating which parental sequence is responsible for each
        alignment column.
    """
  
    n = len(x)  # number of rows
    m = len(y[0])  # number of cols
    S = substitution_matrix
    M = [[[0 for j in range(m+1)] for i in range(n+1)]
         for k in range(2)]  # order: M[k][i][j] for score
    T = [[[0 for j in range(m+1)] for i in range(n+1)]
         for k in range(2)]  # traceback pointer matrix

    # The initialization of first row and col
    for j in range(1, m+1):
        for k in range(2):
            A = (M[k][0][j-1] + S['-', y[k][j-1]], 4)
            B = (M[1-k][0][j-1] + S['-', y[k][j-1]] + r, 1)
            M[k][0][j], T[k][0][j] = max(A, B)

    for i in range(1, n+1):
        for k in range(2):
            A = (M[k][i-1][0] + S[x[i-1], '-'], 5)
            B = (M[1-k][i-1][0] + S[x[i-1], '-'] + r, 2)
            M[k][i][0], T[k][i][0] = max(A, B)

    # Fill matrix table
    for i in range(1, n+1):
        for j in range(1, m+1):
            for k in range(2):
                A = (M[k][i-1][j-1] + S[x[i-1], y[k][j-1]], 6)
                B = (M[k][i-1][j] + S[x[i-1], '-'], 5)
                C = (M[k][i][j-1] + S['-', y[k][j-1]], 4)
                D = (M[1-k][i-1][j-1] + S[x[i-1], y[k][j-1]] + r, 3)
                E = (M[1-k][i-1][j] + S[x[i-1], '-'] + r, 2)
                F = (M[1-k][i][j-1] + S['-', y[k][j-1]] + r, 1)

                M[k][i][j], T[k][i][j] = max(A, B, C, D, E, F)

    score = max(M[0][n][m], M[1][n][m])

    # traceback setup
    alignment = ['', '', '']
    parents = ''

    k = 1
    if M[0][n][m] >= M[1][n][m]:
        k = 0

    i = n
    j = m
    # print(T,k)

    # begin traceback
    while not (i == 0 and j == 0):
        pointer = T[k][i][j]
        parents = str(k) + parents
        if pointer == 6:
            alignment[0] = x[i-1] + alignment[0]
            alignment[1] = y[0][j-1] + alignment[1]
            alignment[2] = y[1][j-1] + alignment[2]
            i -= 1
            j -= 1
        if pointer == 5:
            alignment[0] = x[i-1] + alignment[0]
            alignment[1] = '-' + alignment[1]
            alignment[2] = '-' + alignment[2]
            i -= 1
        if pointer == 4:
            alignment[0] = '-' + alignment[0]
            alignment[1] = y[0][j-1] + alignment[1]
            alignment[2] = y[1][j-1] + alignment[2]
            j -= 1
        if pointer == 3:
            alignment[0] = x[i-1] + alignment[0]
            alignment[1] = y[0][j-1] + alignment[1]
            alignment[2] = y[1][j-1] + alignment[2]
            i -= 1
            j -= 1
            k = 1-k
        if pointer == 2:
            alignment[0] = x[i-1] + alignment[0]
            alignment[1] = '-' + alignment[1]
            alignment[2] = '-' + alignment[2]
            i -= 1
            k = 1-k
        if pointer == 1:
            alignment[0] = '-' + alignment[0]
            alignment[1] = y[0][j-1] + alignment[1]
            alignment[2] = y[1][j-1] + alignment[2]
            j -= 1
            k = 1-k

    #print(score, alignment, parents)
    return (score, alignment, parents)


def pprint_recombine_align(recombine_align_result, names=["x", "y0", "y1"], width=80, stream=sys.stdout):
    """Pretty prints the result of recombine_align.

    Args:
        recombine_align_result: an output from recombine_align, which is a tuple (score, alignment, parents)
        names: a list of three strings giving the names of the three sequences x, y_0, and y_1
        width: an integer specifying the number of columns of the alignment to print before wrapping
        stream: the stream (a file-like object) to which to print
    """
    score, alignment, parents = recombine_align_result
    num_columns = len(parents)
    name_width = max(map(len, names))
    all_names = names + [""]
    all_rows = alignment + [parents]

    print(f"Score: {score}", file=stream)
    for start_column in range(0, num_columns, width):
        for name, row in zip(all_names, all_rows):
            row_slice = row[start_column: start_column + width]
            print(f"{name:>{name_width}} {row_slice}", file=stream)
        print(file=stream)
