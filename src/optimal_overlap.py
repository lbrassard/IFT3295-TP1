import numpy as np
from pathlib import Path
import csv
import sys


MATCH = 4
MISMATCH = -4
INDEL = -8
DIAG = 1
LEFT = 2
UP = 3


def find_optimal_overlap(w: str, v: str):
    # w on top, v on side
    m = len(w)
    n = len(v)

    scores = np.zeros((n+1, m+1), dtype=int)
    paths = np.zeros((n+1, m+1), dtype=int)

    # first num is score, second num is row it's found in
    opt_score = [0, 0]

    # n, i = row;  m, j = column
    for i in range(1, n+1):
        for j in range(1, m+1):
            align = scores[i - 1][j - 1] + MATCH if w[j - 1] == v[i - 1] else scores[i - 1][j - 1] + MISMATCH
            insert = scores[i][j-1] + INDEL
            delete = scores[i-1][j] + INDEL

            # Only keeping one alignment - preference for alignment because it's simpler

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
    alignment = ["", ""]
    num_seq_1 = 0
    num_seq_2 = 0

    while current != 0:
        if current == DIAG:
            alignment[0] = w[current_column - 1] + alignment[0]
            alignment[1] = v[current_row - 1] + alignment[1]
            current_column -= 1
            current_row -= 1
            num_seq_1 += 1
            num_seq_2 += 1
        elif current == LEFT:
            alignment[0] = w[current_column - 1] + alignment[0]
            alignment[1] = "-" + alignment[1]
            current_column -= 1
            num_seq_1 += 1
        else:
            alignment[0] = "-" + alignment[0]
            alignment[1] = v[current_row - 1] + alignment[1]
            current_row -= 1
            num_seq_2 += 1
        current = paths[current_row][current_column]
        length += 1

    rest_seq_1 = len(w) - num_seq_1

    alignment[0] = w[:rest_seq_1] + alignment[0]
    alignment[1] = (" " * rest_seq_1) + alignment[1] + v[num_seq_2:]

    return opt_score[0], alignment, length


def get_two_sequences(file: str):
    sequence1 = ""
    sequence2 = ""

    path = Path(file)
    if not path.is_file():
        print(f"{path} is not a file.")
        return

    with open(path) as file:
        line_num = 1
        sequence1_line = 2
        sequence2_line = 6
        for line in file:
            if line_num == sequence1_line:
                sequence1 += line[:-1]
            if line_num == sequence2_line:
                sequence2 += line[:-1]
            if line_num > 8:
                print("Error: file not in expected format.")
            line_num += 1

    return sequence1, sequence2


def main():
    if len(sys.argv) != 2:
        print("Expected a FASTQ file as a single argument.")
        return

    sequences = get_two_sequences(sys.argv[1])
    results = find_optimal_overlap(sequences[0], sequences[1])
    print("Score: ", end='')
    print(results[0])
    print(results[1][0])
    print(results[1][1])
    print("Length: ", end='')
    print(results[2])


if __name__ == "__main__":
    main()


