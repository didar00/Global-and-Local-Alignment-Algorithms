"""

BBM411: Fundamentals of Bioinformatics (Fall 2021)
Assignment 1
Part 2


"""

import numpy as np

# python main.py --input=input.txt --gap_extension=-5 --gap_opening=-10 --score_matrix=scoring_matrices/BLOSUM62.txt --alg=G

aa_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', "*"]

class Sequence_Alignment:
    
    def __init__(self, P, Q, score_matrix, gap_ex, gap_op, score_mat_file, alignment_alg = 'G', outfile="output.txt"):
        self.P = P
        self.Q = Q
        self.score_matrix = score_matrix
        self.ALGORITHM = alignment_alg
        self.GAP_EX = gap_ex
        self.GAP_OP = gap_op
        self.alignment_matrix = np.zeros((len(Q)+1, len(P)+1))
        self.mapped_amino_acids = {}
        self.outfile = outfile
        self.percent_identity = -333
        self.p_result = ""
        self.q_result = ""
        self.score_mat_file = score_mat_file
        # index each amino acid symbol
        for ind, aa in enumerate(aa_list):
            self.mapped_amino_acids[aa] = ind

        if self.ALGORITHM == "G":
            self.global_alignment()
        elif self.ALGORITHM == "L":
            self.local_alignment()

    def get_percent_identity(self):
        return self.percent_identity


    def global_alignment(self):
        p = self.P
        q = self.Q
        mat = self.alignment_matrix
        trace_mat = [[(-1,1) for j in range(len(p)+1)] for i in range(len(q)+1)]


        # fill in the vertival border values 
        # with gap opening penalty scores
        for i in range(1, len(q)+1):
            if i == 0:
                mat[i][0] = self.GAP_OP
            else:
                mat[i][0] = mat[i-1][0] + self.GAP_EX
            trace_mat[i][0] = (i-1, 0)

        for j in range(1, len(p)+1):
            if j == 0:
                mat[0][j] = self.GAP_OP
            else:
                mat[0][j] = mat[0][j-1] + self.GAP_EX
            trace_mat[0][j] = (0,j-1)


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

        while (i,j) != (0,0):
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

        # number of matches in alignment
        match_count = 0
        # alignment representation
        # to mark the alignment between
        # sequences
        align_rep = []
        for i, j in zip(p_result, q_result):
            if i == j: # match
                match_count += 1
            if i != j: # no match
                align_rep.append(" ")
            else: # match
                align_rep.append("|")


        self.percent_identity = match_count * 100 / len(p_result)


        data = ["Global Alignment Results (Needleman-Wunch)\n",
                    "Gap opening penalty  :{} \n".format(self.GAP_OP),
                    "Gap extension penalty  :{} \n".format(self.GAP_EX),
                    "Score matrix   :{} \n".format(self.score_mat_file),
                    "".join(p_result), "\n", 
                    "".join(align_rep), "\n",
                    "".join(q_result), "\n",
                    "Raw alignment score  :{} \n".format(mat[-1][-1]),
                    "Match rate :{} percent \n".format(self.percent_identity),"\n\n"]
        self.write_output_to_file(data)

        self.p_result = p_result
        self.q_result = q_result
        return p_result, q_result


    def local_alignment(self):
        p = self.P
        q = self.Q
        mat = self.alignment_matrix
        trace_mat = [[(-1,-1) for j in range(len(p)+1)] for i in range(len(q)+1)]
        trace_mat[0][0] = (0,0)

        # fill in the vertival border values 
        # with gap opening penalty scores
        for i in range(1, len(q)+1):
            if i == 0:
                mat[i][0] = self.GAP_OP
            else:
                mat[i][0] = mat[i-1][0] + self.GAP_EX
            trace_mat[i][0] = (i-1, 0)

        for j in range(1, len(q)+1):
            if j == 0:
                mat[0][j] = self.GAP_OP
            else:
                mat[0][j] = mat[0][j-1] + self.GAP_EX
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
            elif mat[x][y] == 0:
                break
            (x, y) = trace_mat[x][y]
        
        # number of matches in alignment
        match_count = 0
        # alignment representation
        # to mark the alignment between
        # sequences
        align_rep = []

        for i, j in zip(p_result, q_result):
            if i == j: # match
                match_count += 1
            if i != j: # no match
                align_rep.append(" ")
            else: # match
                align_rep.append("|")
        
        
        self.percent_identity = match_count * 100 / len(p_result)
        data = ["Local Alignment Results (Smith-Waterman)\n",
                    "Gap opening penalty  :{} \n".format(self.GAP_OP),
                    "Gap extension penalty  :{} \n".format(self.GAP_EX),
                    "Score matrix   :{} \n".format(self.score_mat_file),
                    "".join(p_result), "\n", 
                    "".join(align_rep), "\n",
                    "".join(q_result), "\n",
                    "Raw alignment score  :{} \n".format(max_),
                    "Match rate :{} percent \n\n".format(self.percent_identity)]
        self.write_output_to_file(data)
        self.p_result = p_result
        self.q_result = q_result
        
        return p_result, q_result

    def write_output_to_file(self, data):
        with open(self.outfile, "a") as f:
            for line in data:
                f.write(line)



    def get_aa_index(self, p_aa, q_aa):
        """
        Returns: the index of amino acid representation letters
        Replace aa's name with * if the amino acid's name is not in the list
        """
        if p_aa not in aa_list:
            p_aa = '*'
        if q_aa not in aa_list:
            q_aa = '*'
        return self.mapped_amino_acids[p_aa], self.mapped_amino_acids[q_aa]



def get_sequences(filename):
    """
    Reads biomoleculer sequences from file
    and returns them
    """
    with open(filename, "r") as f:
        lines = f.readlines()

    P, Q = lines[0].strip(), lines[1].strip()
    return P, Q



def get_score_matrix(mat_file):
    """
    Read the score matrix from file
    and store it in numpy array
    """
    with open(mat_file, "r") as f:
        lines = [line.rstrip() for line in f]

    n = len(lines)-1
    
    score_matrix = np.zeros((n,n))
    for i in range(1, n+1):
        line = lines[i].strip().split(" ")[1:]
        line = [ele for ele in line if ele.strip()]
        for j in range(n):
            score_matrix[i-1][j] = int(line[j])
    return score_matrix
            

def input_handler(args):
    """
    Prepares the input arguments for sequence alignment
    """
    score_matrix = get_score_matrix(args.score_matrix)

    P, Q = get_sequences(args.input)
    SA = Sequence_Alignment(P, Q, score_matrix, int(args.gap_extension), int(args.gap_opening), args.score_matrix, args.alg)


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
    parser.add_argument('--alg', required=False,
                        help='Alignment algorithm: L for local; G for global (default)')
    parser.add_argument('--output_file_name', required=False,
                        metavar="/path/to/output/file/",
                        help="File name of the output file")
    args = parser.parse_args()
    input_handler(args)


