import numpy as np

# O, U, J deleted
aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X']

class Sequence_Alignment:
    
    def __init__(self, P, Q, score_matrix, alignment_alg, gap_ex, gap_op):
        self.P = P
        self.Q = Q
        self.score_matrix = score_matrix
        self.ALGORITHM = alignment_alg
        self.GAP_EX = gap_ex
        self.GAP_OP = gap_op
        self.alignment_matrix = np.zeros((len(Q)+1, len(P)+1))
        self.mapped_amino_acids = {}
        # index each amino acid symbol
        for ind, aa in enumerate(aa_list):
            self.mapped_amino_acids[aa] = ind

        if self.ALGORITHM == "G":
            self.global_alignment()
        elif self.ALGORITHM == "L":
            self.local_alignment()

    def global_alignment(self):
        p = self.P
        q = self.Q
        print(p, q)
        mat = self.alignment_matrix
        #backtracking_mat = np.zeros((len(p), len(q)))
        trace_mat = [[(-1,1) for j in range(len(p)+1)] for i in range(len(q)+1)]
        print(trace_mat)


        # fill in the vertival border values 
        # with gap opening penalty scores
        for i in range(1, len(q)+1):
            if i == 0:
                mat[i][0] = self.GAP_OP
            else:
                mat[i][0] = mat[i-1][0] + self.GAP_EX # NOT SURE IF THAT WHAT IT SHOULD BE!!!!!
            trace_mat[i][0] = (i-1, 0)

        for j in range(1, len(p)+1):
            if j == 0:
                mat[0][j] = self.GAP_OP
            else:
                mat[0][j] = mat[0][j-1] + self.GAP_EX
            trace_mat[0][j] = (0,j-1)
        print(mat)
        for line in trace_mat:
            print(line)

        # start to fill in the sequence matrix
        # row by row
        for i in range(1, len(q)+1):
            for j in range(1, len(p)+1):
                # get the indices of related amino acids
                # to obtain the scores of them in the score matrix
                p_ind, q_ind = self.get_aa_index(p[j-1], q[i-1])
                # get the maximum value out of 3 possible neighbors


                max_ = max(mat[i-1][j] + self.GAP_OP, mat[i][j-1] + self.GAP_OP, mat[i-1][j-1] + self.score_matrix[p_ind][q_ind])
                mat[i][j] = max_
                #print("***", i , j)

                # backtracking, but a scarce matrix, too much space!!!!
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!12
                if mat[i][j] == mat[i-1][j] + self.GAP_OP:
                    trace_mat[i][j] = (i-1, j)
                    # check if the previous element was also a gap
                    if i > 1 and trace_mat[i-1][j] == (i-2, j):
                        # replace gap opening penalty with gap extension penalty
                        mat[i][j] = mat[i][j] - self.GAP_OP + self.GAP_EX
                elif mat[i][j] == mat[i][j-1] + self.GAP_OP:
                    trace_mat[i][j] = (i, j-1)
                    # check if the previous element was also a gap
                    if j > 1 and trace_mat[i][j-1] == (i, j-2):
                        # replace gap opening penalty with gap extension penalty
                        mat[i][j] = mat[i][j] - self.GAP_OP + self.GAP_EX
                elif mat[i][j] == mat[i-1][j-1] + self.score_matrix[p_ind][q_ind]:
                    trace_mat[i][j] = (i-1, j-1)



        # resulting alignments of the sequences
        p_result = ""
        q_result = ""

        # start backtracking
        (i, j) = (len(q), len(p))
        print(mat)

        for line in trace_mat:
            print(line)

        #k = 0
        while (i,j) != (0,0):
            #print(p_result, q_result, i, j)
            if trace_mat[i][j] == (i-1, j-1):
                p_result =  p[j-1] + p_result
                q_result = q[i-1] + q_result
            elif trace_mat[i][j] == (i-1, j):
                p_result = "-" + p_result
                q_result = q[i-1] + q_result
            elif trace_mat[i][j] == (i, j-1):
                p_result = p[j-1] + p_result
                q_result = "-" + q_result
            (i, j) = trace_mat[i][j]
            """ if k < 10:
                print(p_result, q_result, i, j)
                k += 1 """


        print(p_result, q_result)
        return p_result, q_result


    def local_alignment(self):
        p = self.P
        q = self.Q
        print(p, q)
        mat = self.alignment_matrix
        trace_mat = [[(-1,-1) for j in range(len(p)+1)] for i in range(len(q)+1)]
        #mat[0][0] = 0
        trace_mat[0][0] = (0,0)

        # set border values of the trace matrix
        for i in range(1, len(q)+1):
            #mat[i][0] = 0
            trace_mat[i][0] = (i-1, 0)

        for j in range(1, len(q)+1):
            #mat[0][j] = 0
            trace_mat[0][j] = (0, j-1)

        for i in range(1, len(q)+1):
            for j in range(1, len(p)+1):
                # get the indices of the related amino acids' letters
                # to obtain scores from the score matrix
                p_ind, q_ind = self.get_aa_index(p[j-1], q[i-1])

                # get the maximum values out of 3 neighbors
                max_ = max(0, mat[i-1][j] + self.GAP_OP, mat[i][j-1] + self.GAP_OP, mat[i-1][j-1] + self.score_matrix[p_ind][q_ind])
                mat[i][j] = max_

                # update the trace matrix
                if mat[i][j] == mat[i-1][j] + self.GAP_OP:
                    trace_mat[i][j] = (i-1, j)
                    # check if the previous element was also a gap
                    if i > 1 and trace_mat[i-1][j] == (i-2, j):
                        # replace gap opening penalty with gap extension penalty
                        mat[i][j] = mat[i][j] - self.GAP_OP + self.GAP_EX
                elif mat[i][j] == mat[i][j-1] + self.GAP_OP:
                    trace_mat[i][j] = (i, j-1)
                    # check if the previous element was also a gap
                    if j > 1 and trace_mat[i][j-1] == (i, j-2):
                        # replace gap opening penalty with gap extension penalty
                        mat[i][j] = mat[i][j] - self.GAP_OP + self.GAP_EX
                elif mat[i][j] == mat[i-1][j-1] + self.score_matrix[p_ind][q_ind]:
                    trace_mat[i][j] = (i-1, j-1)

         
        
        print(mat)

        for line in trace_mat:
            print(line)


        """
        find the maximum score and coordinates
        in alignment matrix starting from bottom-right
        """

        max_ = 0 # minimum value can be 0 in local alignment

        for i in range(len(q)+1):
            for j in range(len(p)+1):
                if mat[i][j] > max_:
                    max_ = mat[i][j]
                    x, y = i, j

        """
        Start backtracing from the coordinates that 
        is found above using maximum aligment score
        """

        # resulting alignments of the sequences
        p_result = ""
        q_result = ""

        while (x,y) != (0,0):
            if trace_mat[x][y] == (x-1, y-1):
                p_result =  p[y-1] + p_result
                q_result = q[x-1] + q_result
            elif trace_mat[x][y] == (x-1, y):
                p_result = "-" + p_result
                q_result = q[x-1] + q_result
            elif trace_mat[x][y] == (x, y-1):
                p_result = p[y-1] + p_result
                q_result = "-" + q_result
            print(x,y)
            (x, y) = trace_mat[x][y]

        print(p_result, q_result)
        return p_result, q_result



    def get_aa_index(self, p_aa, q_aa):
        """
        Returns: the index of amino acid representation letters
        """
        return self.mapped_amino_acids[p_aa], self.mapped_amino_acids[q_aa]
        
        



def get_sequences(filename):
    """
    Reads biomoleculer sequences from file
    and returns them
    """
    with open(filename, "r") as f:
        lines = f.readlines()
    #print(lines)
    P, Q = check_seq(lines[0].strip()), check_seq(lines[1].strip())

    return P, Q


def check_seq(seq):
    """
    Checks whether the amino acid sequence is valid or not
    """
    for aa in seq:
        if aa not in aa_list:
            raise Exception("Invalid amino acid!")
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
    """
    Prepares the input arguments for sequence alignment
    """
    P, Q = get_sequences(args.input)
    score_matrix = get_score_matrix(args.score_matrix)
    SA = Sequence_Alignment(P, Q, score_matrix, args.alg, int(args.gap_extension), int(args.gap_opening))
    
    




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


