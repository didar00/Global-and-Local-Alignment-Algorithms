import numpy as np

aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'O', 'S', 'U', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', 'J']

class Sequence_Alignment:
    
    def __init__(self, P, Q, score_matrix, alignment_alg, gap_ex, gap_pen):
        self.P = P
        self.Q = Q
        self.score_matrix = score_matrix
        self.ALGORITHM = alignment_alg
        self.GAP_EX = gap_ex
        self.GAP_PEN = gap_pen




def get_sequences(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    P, Q = check_seq(lines[0]), check_seq(lines[1])


    return P, Q

def check_seq(seq):
    for aa in seq:
        if aa not in aa_list:
            return False
    return seq
            

def input_handler(args):
    P, Q = get_sequences(args.input)
    SA = Sequence_Alignment(P, Q, args.score_matrix, args.alg, args.gap_extension, args.gap_penalty)




if __name__ == '__main__':
    import argparse

    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Sequence alignment')
    parser.add_argument('--input', required=True,
                        metavar="/path/to/input/file/",
                        help='File name/path of the input file')
    parser.add_argument('--score_matrix', required=True,
                        metavar="path or URL to score_matrix",
                        help='File name for score matrix')
    parser.add_argument('--gap_extension', required=True,
                        help='Gap extention score')
    parser.add_argument('--gap_penalty', required=True,
                        help='Gap penalty score')
    parser.add_argument('--alg', required=True,
                        help='Alignment algorithm: L for local; G for global')

    args = parser.parse_args()
    input_handler(args)


