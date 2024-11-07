"""
Short read aligner based on Burrows–Wheeler transform.

Author :
    Gloria BENOIT
Date :
    2024-10-11
"""

import sys
import time

def get_BWT(seq):
    """ Get Burrows–Wheeler transform.

    Parameters
    ----------
    seq : str
        Sequence on which you wish to align reads.

    Returns
    -------
    str
        Burrows–Wheeler transform.
    """
    new_seq = '$' + seq
    circ_perm = []

    # Construction de la matrice
    for i in range(len(new_seq)):
        circ_perm.append(new_seq)
        new_seq = new_seq[-1] + new_seq[:-1]

    # Ordonnée de la matrice
    tri = sorted(circ_perm)
    bwt = ''
    for ligne in tri:
        bwt += ligne[-1]

    return bwt

def get_indexes(bwt):
    """ Get position indexes and cumulative indexes.

    Parameters
    ----------
    bwt : str
        Burrows–Wheeler transform.

    Returns
    -------
    list
        Position indexes.
    dict
        Cumulative indexes, with nucleotide as keys and starting position
        index as values.
    """
    all_values = sorted(list(set(bwt)))
    num_id = {}
    cumul = {}

    for value in all_values:
        num_id[value] = 0
        cumul[value] = 0

    num = []
    for char in bwt:
        num.append(num_id[char])
        num_id[char] += 1

    prev_count = 0
    for char, count in num_id.items():
        cumul[char] = prev_count
        prev_count += count

    return num, cumul

def get_origin_pos(bwt):
    """ Get position in original sequence.

    Parameters
    ----------
    bwt : str
        Burrows–Wheeler transform.

    Returns
    -------
    list
        Position indexes.
    """
    num, cumul = get_indexes(bwt)
    L = len(bwt)
    pos_origin = [0] * L
    pos = L - 1
    i = 0

    while pos >= 0:
        pos_origin[i] = pos
        i = cumul[bwt[i]] + num[i]
        pos = pos - 1

    return pos_origin

def get_FM(bwt):
    """ Get FM-index.

    Parameters
    ----------
    bwt : str
        Burrows–Wheeler transform.

    Returns
    -------
    matrix
        FM-index.
    """
    FM = []
    all_values = sorted(list(set(bwt)))
    current_dict = {}
    for value in all_values:
        current_dict[value] = 0
    
    for char in bwt:
        current_dict = current_dict.copy()
        current_dict[char] += 1
        FM.append(current_dict)
    
    FM.append(dict(zip(all_values, [0] * len(bwt))))
    return FM

def mapping(sequence, read, bwt):
    """ Search for read in sequence.

    Parameters
    ----------
    sequence : str
        Sequence in which you wish to search.
    read : str
        Sequence you wish to search.
    bwt : str
        Burrows–Wheeler transform.

    Returns
    -------
    int
        Beginning of position.
    int
        End of position.
    """
    num, cumul = get_indexes(bwt)
    FM = get_FM(bwt)

    i = len(read) - 1
    nt = read[i]

    b = cumul[nt]
    e = cumul[nt] + FM[len(bwt) - 1][nt] - 1
    
    while (i > 0 and b <= e):
        nt = read[i - 1]
        b = cumul[nt] + FM[b - 1][nt]
        e = cumul[nt] + FM[e][nt] - 1
        i = i - 1
    if e >= b:
        return b, e
    else:
        sys.exit(f"{read} not found in sequence.")

def get_reads_pos(begin, end, pos_origin):
    """ Get reads positions.

    Parameters
    ----------
    begin : int
        Beginning of position.
    end : int
        End of position.
    pos_origin : list
        Position indexes.

    Returns
    -------
    list
        Position of found reads in original sequence.
    """
    reads = []
    for i in range(begin, end + 1):
        reads.append(pos_origin[i] + 1)
    return reads

def read_fasta(file):
    """ Read fasta file.

    Parameters
    ----------
    file
        Fasta file.

    Returns
    -------
    str
        Sequence.
    """
    with open(file, 'r') as f_in:
        sequence = ""
        for ligne in f_in:
            if ligne.startswith(">"):
                continue
            else:
                sequence += ligne.strip()

        return sequence

def mapper(sequence, read):
    """ Find reads in sequence.

    Parameters
    ----------
    sequence : str
        Sequence in which you wish to search.
    read : str
        Sequence you wish to search.
    """
    bwt = get_BWT(sequence)
    pos_origin = get_origin_pos(bwt)
    
    begin, end = mapping(sequence, read, bwt)
    reads = get_reads_pos(begin, end, pos_origin)
    print(f"Found {len(reads)} matchs:")
    for read in sorted(reads):
        print(read)


if __name__ == "__main__":
    # Récupération des arguments
    if len(sys.argv) != 3:
        sys.exit("Wrong number of arguments!\nUsage: mapper.py sequence.fasta short_read")

    SEQUENCE = read_fasta(sys.argv[1])
    READ = sys.argv[2]

    # Mapping
    start_time = time.time()
    mapper(SEQUENCE, READ)
    print(f"Execution: {(time.time() - start_time)} seconds")