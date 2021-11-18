import numpy as np

aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'O', 'S', 'U', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', 'J']

class Sequence_Alignment:
    
    def __init__(self, P, Q, score_matrix, alignment_alg, gap_ex, gap_op):
        self.P = P
        self.Q = Q
        self.score_matrix = score_matrix
        self.ALGORITHM = alignment_alg
        self.GAP_EX = gap_ex
        self.GAP_OP = gap_op
        self.matrix = np.zeros((len(Q)+1, len(P)+1))
        self.mapped_amino_acids = {}
        for ind, aa in enumerate(aa_list):
            #print(ind, aa)
            self.mapped_amino_acids[aa] = ind

    def global_alignment(self):
        p = self.P
        q = self.Q
        print(p, q)
        mat = self.matrix
        #backtracking_mat = np.zeros((len(p), len(q)))
        backtracking_mat = [[(-1,1) for j in range(len(p)+1)] for i in range(len(q)+1)]
        print(backtracking_mat)


        # fill in the vertival border values 
        # with gap opening penalty scores
        for i in range(1, len(q)+1):
            if i == 0:
                mat[i][0] = self.GAP_OP
            else:
                mat[i][0] = mat[i-1][0] + self.GAP_EX # NOT SURE IF THAT WHAT IT SHOULD BE!!!!!
            backtracking_mat[i][0] = (i-1, 0)

        for j in range(1, len(p)+1):
            if j == 0:
                mat[0][j] = self.GAP_OP
            else:
                mat[0][j] = mat[0][j-1] + self.GAP_EX
            backtracking_mat[0][j] = (0,j-1)
        print(mat)
        print(backtracking_mat)

        # start to fill in the sequence matrix
        # row by row
        for i in range(1, len(q)+1):
            for j in range(1, len(p)+1):
                # get the indices of related amino acids
                # to obtain the scores of them in the score matrix
                p_ind, q_ind = self.get_aa_index(p[j-1], q[i-1])
                # get the maximum value out of 3 possible neighbors


                max_ = max(mat[i-1][j] + self.GAP_OP, mat[i][j-1] + self.GAP_OP, mat[i-1][j-1] + self.score_matrix[p_ind][q_ind])
                print(max_)
                mat[i][j] = max_
                #print("***", i , j)

                # backtracking, but a scarce matrix, too much space!!!!
                if mat[i][j] == mat[i-1][j] + self.GAP_OP:
                    backtracking_mat[i][j] = (i-1, j)
                    # check if the previous element was also a gap
                    if i > 1 and backtracking_mat[i-1][j] == (i-2, j):
                        # replace gap opening penalty with gap extension penalty
                        mat[i][j] = mat[i][j] - self.GAP_OP + self.GAP_EX
                elif mat[i][j] == mat[i][j-1] + self.GAP_OP:
                    backtracking_mat[i][j] = (i, j-1)
                    # check if the previous element was also a gap
                    if j > 1 and backtracking_mat[i][j-1] == (i, j-2):
                        # replace gap opening penalty with gap extension penalty
                        mat[i][j] = mat[i][j] - self.GAP_OP + self.GAP_EX
                elif mat[i][j] == mat[i-1][j-1] + self.score_matrix[p_ind][q_ind]:
                    backtracking_mat[i][j] = (i-1, j-1)



        # resulting alignments of the sequences
        p_result = ""
        q_result = ""

        # start backtracking
        (i, j) = (len(q), len(p))
        print(mat)

        for line in backtracking_mat:
            print(line)

        k = 0
        while (i,j) != (0,0):
            #print(p_result, q_result, i, j)
            if backtracking_mat[i][j] == (i-1, j-1):
                p_result =  p[j-1] + p_result
                q_result = q[i-1] + q_result
            elif backtracking_mat[i][j] == (i-1, j):
                p_result = "-" + p_result
                q_result = q[i-1] + q_result
            elif backtracking_mat[i][j] == (i, j-1):
                p_result = p[j-1] + p_result
                q_result = "-" + q_result
            (i, j) = backtracking_mat[i][j]
            if k < 10:
                print(p_result, q_result, i, j)
                k += 1

            #print(i,j)

        print(p_result, q_result)
        return p_result, q_result


    def local_alignment(self):
        pass

    def get_aa_index(self, p_aa, q_aa):
        """
        Returns: the index of amino acid representation letters
        """
        return self.mapped_amino_acids[p_aa], self.mapped_amino_acids[q_aa]
        
        





def get_sequences(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    #print(lines)
    P, Q = check_seq(lines[0].strip()), check_seq(lines[1].strip())


    return P, Q


def check_seq(seq):
    for aa in seq:
        if aa not in aa_list:
            return None
    return seq


def get_score_matrix(mat_file):
    """
    Read the score matrix from file
    and store it in numpy array
    """
    with open(mat_file, "r") as f:
        lines = [line.rstrip() for line in f]
        #lines = f.readlines()

    n = len(lines)-1
    #print(lines)
    
    score_matrix = np.zeros((n,n))
    for i in range(1, n+1):
        line = lines[i].strip().split(" ")[1:]
        line = [ele for ele in line if ele.strip()]
        #print(line)
        for j in range(n):
            score_matrix[i-1][j] = int(line[j])

    return score_matrix
            

def input_handler(args):
    P, Q = get_sequences(args.input)
    score_matrix = get_score_matrix(args.score_matrix)
    SA = Sequence_Alignment(P, Q, score_matrix, args.alg, int(args.gap_extension), int(args.gap_opening))
    SA.global_alignment()
    




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
                        help='Gap extention penalty score')
    parser.add_argument('--gap_opening', required=True,
                        help='Gap opening penalty score')
    parser.add_argument('--alg', required=True,
                        help='Alignment algorithm: L for local; G for global')

    args = parser.parse_args()
    input_handler(args)


