import numpy as np

# python main.py --input=input.txt --gap_extension=-1 --gap_opening=-1 --score_matrix=blosum62.txt --alg=G

# O, U, J deleted
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
        #print(p, q)
        mat = self.alignment_matrix
        #backtracking_mat = np.zeros((len(p), len(q)))
        trace_mat = [[(-1,1) for j in range(len(p)+1)] for i in range(len(q)+1)]
        #print(trace_mat)


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
        """ print(mat)
        for line in trace_mat:
            print(line) """

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
        """ print(mat)

        for line in trace_mat:
            print(line) """

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
        """ print("*************************")
        print(mat)
        print(p_result)
        print(q_result)
        print(match_count)
        print(len(p_result))
        print("*************************") """

        self.percent_identity = match_count * 100 / len(p_result)


        data = ["Global Alignment Results (Needleman-Wunch)\n",
                    "Gap opening penalty  :{} \n".format(self.GAP_OP),
                    "Gap extension penalty  :{} \n".format(self.GAP_EX),
                    "Score matrix   :{} \n".format(self.score_mat_file),
                    "".join(p_result), "\n", 
                    "".join(align_rep), "\n",
                    "".join(q_result), "\n",
                    "Raw alignment score  :{} \n".format(mat[-1][-1]),
                    "Match rate :{} percent \n".format(self.percent_identity),"\n"]
        self.write_output_to_file(data)

        self.p_result = p_result
        self.q_result = q_result
        return p_result, q_result


    def local_alignment(self):
        p = self.P
        q = self.Q
        #print(p, q)
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


                """
                What if none of them is zero but there is a zero??
                ERROR!!!!!!!!!!!!!!!!!!!!!!!!
                
                """
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

         
        
        """ print(mat)

        for line in trace_mat:
            print(line) """


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
            #print(x,y)
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
        

        print(p_result)
        print(q_result)
        print(match_count)
        
        self.percent_identity = match_count * 100 / len(p_result)
        data = ["Local Alignment Results (Smith-Waterman)\n",
                    "Gap opening penalty  :{} \n".format(self.GAP_OP),
                    "Gap extension penalty  :{} \n".format(self.GAP_EX),
                    "Score matrix   :{} \n".format(self.score_mat_file),
                    "".join(p_result), "\n", 
                    "".join(align_rep), "\n",
                    "".join(q_result), "\n",
                    "Raw alignment score  :{} \n".format(max_),
                    "Match rate :{} percent \n".format(self.percent_identity)]
        self.write_output_to_file(data)
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
        
        


class MSA:

    def __init__(self, sequence_list, score_matrix, gap_ex, gap_op, outfile="output.txt"):
        self.sequence_list = sequence_list
        self.seq_count = len(sequence_list)
        #print("SEQ COUNT", self.seq_count)
        print("sequence list" , sequence_list)
        #self.similarity_mat = np.zeros((self.seq_count, self.seq_count))
        self.similarity_matrix = dict()
        self.score_matrix = score_matrix
        self.GAP_EX = gap_ex
        self.GAP_OP = gap_op
        self.outfile = outfile
        self.aligned_seqs = list()

        # construct the similarity matrix
        self.fill_similarity_matrix()
        # align multiple sequences
        self.align_seqs()
        # print the aligned sequences one after another
        self.print_MSA()


    def align_seqs(self):

        aligned_seqs_names = list()
        aa = dict()

        for p, q in self.similarity_matrix:
            if p in aligned_seqs_names and q in aligned_seqs_names:
                continue
            elif p in aligned_seqs_names:
                aligned_seqs_names.append(p)
                aligned_seqs_names.append(q)
                p = aa[p]
                SA = Sequence_Alignment(p, q, self.score_matrix, self.GAP_EX, self.GAP_OP)
                print("resulting p,q ", SA.p_result, SA.q_result)
                consensus_seq = self.create_consensus(SA.p_result, SA.q_result)
                #self.similarity_matrix[(consensus_seq,q)] = self.similarity_matrix[(p,q)]
                #del self.similarity_matrix[(p,q)]
                self.aligned_seqs.append(SA.q_result)
                aa[p] = consensus_seq
                aa[q] = consensus_seq
            elif q in aligned_seqs_names:
                aligned_seqs_names.append(p)
                aligned_seqs_names.append(q)
                q = aa[q]
                SA = Sequence_Alignment(p, q, self.score_matrix, self.GAP_EX, self.GAP_OP)
                print("resulting p,q ", SA.p_result, SA.q_result)
                consensus_seq = self.create_consensus(SA.p_result, SA.q_result)
                self.aligned_seqs.append(SA.p_result)
                aa[p] = consensus_seq
                aa[q] = consensus_seq
                #self.similarity_matrix[(p,consensus_seq)] = self.similarity_matrix[(p,q)]
                #del self.similarity_matrix[(p,q)]
            else:
                aligned_seqs_names.append(p)
                aligned_seqs_names.append(q)
                SA = Sequence_Alignment(p, q, self.score_matrix, self.GAP_EX, self.GAP_OP)
                print("resulting p,q ", SA.p_result, SA.q_result)
                consensus_seq = self.create_consensus(SA.p_result, SA.q_result)
                self.aligned_seqs.append(SA.p_result)
                self.aligned_seqs.append(SA.q_result)
                aa[p] = consensus_seq
                aa[q] = consensus_seq
            #print(p,q)
            print("consensus ==>", consensus_seq)
            
            """
            print(SA.p_result)
            print(SA.q_result)
            """

            #self.aligned_seqs.append(consensus_seq)

            aligned_seqs_names.append(p)
            aligned_seqs_names.append(q)




    def fill_similarity_matrix(self):
        """
        Constructs the similarity matrix using global/local
        pairwise alignment
        Returns : similarity matrix
        """
        for i in range(1, self.seq_count):
            for j in range(0, i):
                if i != j:
                    P = self.sequence_list[i]
                    Q = self.sequence_list[j]
                    SA = Sequence_Alignment(P, Q, self.score_matrix, int(args.gap_extension), int(args.gap_opening), "G")
                    self.similarity_matrix[(P,Q)] = SA.get_percent_identity()
        
        self.similarity_matrix = dict(sorted(self.similarity_matrix.items(), key=lambda item: item[1]))
        print("SIMILARITY MATRIX" , self.similarity_matrix)

    def print_MSA(self):
        print("BEGIN")
        for line in self.aligned_seqs:
            print(line)
        print("END")

    
    def create_consensus(self, seq1, seq2):
        """Sequences must be strings, have the same length, and be aligned"""
        out_seq = ""
        for i, nucleotide in enumerate(seq1):
            couple = [nucleotide, seq2[i]]
            if couple[0] == "-":
                out_seq += couple[1]
            elif couple[1] == "-":
                out_seq += couple[0]
            elif couple[0] == couple[1]:
                out_seq += couple[0]
            elif not couple[0] == couple[1]:
                out_seq += couple[0]
        return out_seq








def get_sequences(filename, task):
    """
    Reads biomoleculer sequences from file
    and returns them
    """
    with open(filename, "r") as f:
        lines = f.readlines()
    #print(lines)
    if task == 1:
        P, Q = lines[0].strip(), lines[1].strip()
        return P, Q
    elif task == 2:
        seq_list = list()
        for line in lines:
            seq_list.append(line.strip())
        return seq_list

    return None


""" def check_seq(seq):
    
    Checks whether the amino acid sequence is valid or not
   
    for aa in seq:
        if aa not in aa_list:
            raise Exception("Invalid amino acid!")
    return seq
 """

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
        
    """ print("****************************************************************************")
    print(score_matrix)
    print("****************************************************************************") """
    return score_matrix
            

def input_handler(args):
    """
    Prepares the input arguments for sequence alignment
    """
    task = int(args.task)
    score_matrix = get_score_matrix(args.score_matrix)
    if task == 1:
        P, Q = get_sequences(args.input, task)
        SA = Sequence_Alignment(P, Q, score_matrix, int(args.gap_extension), int(args.gap_opening), args.score_matrix, args.alg)
    elif task == 2:
        seq_list = get_sequences(args.input, task)
        MSA_ = MSA(seq_list, score_matrix, int(args.gap_extension), int(args.gap_opening))
    



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
                        help='Alignment algorithm: L for local; G for global')
    parser.add_argument('--output_file_name', required=False,
                        metavar="/path/to/output/file/",
                        help="File name of the output file")
    parser.add_argument('--task', required=True,
                        help="Task number: 1 for global_local alignment, 2 for ClustalW")

    args = parser.parse_args()
    input_handler(args)


