import numpy as np

MATCH = 4
MISMATCH = -4
INDEL = -8
DIAG = 1
LEFT = 2
UP = 3


def optimal_overlap(w: str, v: str):
    # w on top, v on side
    m = len(w)
    n = len(v)

    scores = np.zeros((n+1, m+1), dtype=int)
    paths = np.zeros((n+1, m+1), dtype=int)

    # first num is score, second num is row it's found in
    opt_score = [0, 0]

    # n, i = row;  m, j = column
    for j in range(1, m+1):
        for i in range(1, n+1):
            align = scores[i - 1][j - 1] + MATCH if w[j - 1] == v[i - 1] else scores[i - 1][j - 1] + MISMATCH
            insert = scores[i][j-1] + INDEL
            delete = scores[i-1][j] + INDEL

            # Only keeping one alignment - preference for alignment because it's easier :P

            options = [align, insert, delete]
            max_value = max(options)
            max_index = options.index(max_value)
            if max_index == 0:
                scores[i][j] = align
                paths[i][j] = DIAG
            elif max_index == 1:
                scores[i][j] = insert
                paths[i][j] = LEFT
            else:
                scores[i][j] = delete
                paths[i][j] = UP

            # last column - update optimal score so far
            if j == m:
                # make sure opt_score never stays at [0, 0] if all scores are below 0
                if i == 1:
                    opt_score[0] = max_value
                    opt_score[1] = i
                elif max_value > opt_score[0]:
                    opt_score[0] = max_value
                    opt_score[1] = i

    # deduce overlap
    current_row = opt_score[1]
    current_column = m
    current = paths[current_row][current_column]
    length = 0
    # w is first, v is second
    overlap = ["", ""]

    while current != 0:
        if current == DIAG:
            overlap[0] = w[current_column - 1] + overlap[0]
            overlap[1] = v[current_row - 1] + overlap[1]
            current_column -= 1
            current_row -= 1
        elif current == LEFT:
            overlap[0] = w[current_column - 1] + overlap[0]
            overlap[1] = "-" + overlap[1]
            current_column -= 1
        else:
            overlap[0] = "-" + overlap[0]
            overlap[1] = v[current_row - 1] + overlap[1]
            current_row -= 1
        current = paths[current_row][current_column]
        length += 1

    # overlap - just the overlap, not the full words for now
    return [[opt_score[0]], [overlap], [length]]


w1 = "CATCCTTCT"
w2 = "CCTTTCACC"
w3 = "CGAATTCGG"
w4 = "ATCGTTGGT"
print(optimal_overlap(w1, w2))
print(optimal_overlap(w3, w4))
