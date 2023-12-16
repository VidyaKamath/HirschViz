from src.score_functions import get_dna_score_func
import numpy as np
import time
import matplotlib.pyplot as plt


UP = (-1,0)
LEFT = (0, -1)
TOPLEFT = (-1, -1)
ORIGIN = (0, 0)

class NeedlemanWunsch(object):
    def __init__(self,v,w):
        self.v = v
        self.w = w
        self.score_func = get_dna_score_func()
        
        # Output
        self.alignment = []
        self.global_score = None
        self.dp_table = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
        self.pointers = [[ORIGIN for j in range(len(w)+1)] for i in range(len(v)+1)]
        self.reported_cells = []

        # time
        self.alignment_time = None
        self.score_time = None
    
    @staticmethod
    def compute_score_quad_space(v, w, score_func):
        m = len(v)
        n = len(w)
        dp_table = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
        for i in range(0, m+1):
            for j in range (0, n+1):
                if(i==0 and j==0):
                    dp_table[i][j] = 0
                elif(i==0 and j>0):
                    dp_table[i][j] = dp_table[i][j-1] + score_func["-"][w[j-1]]
                elif(i>0 and j==0):
                    dp_table[i][j] = dp_table[i-1][j] + score_func[v[i-1]]["-"]
                else:
                    match_or_mismatch = dp_table[i-1][j-1] + score_func[v[i-1]][w[j-1]]
                    deletion = dp_table[i-1][j] + score_func[v[i-1]]["-"]
                    insertion = dp_table[i][j-1] + score_func["-"][w[j-1]]
                    dp_table[i][j] = max(match_or_mismatch, deletion, insertion)
        return dp_table[m][n]
    
    def __traceback_global(self):
        i,j = len(self.v), len(self.w)
        new_v = []
        new_w = []
        self.reported_cells.append((i,j))
        while True:
            di, dj = self.pointers[i][j]
            if (di,dj) == LEFT:
                new_v.append('-')
                new_w.append(self.w[j-1])
            elif (di,dj) == UP:
                new_v.append(self.v[i-1])
                new_w.append('-')
            elif (di,dj) == TOPLEFT:
                new_v.append(self.v[i-1])
                new_w.append(self.w[j-1])
            i, j = i + di, j + dj
            self.reported_cells.append((i, j))
            if (i <= 0 and j <= 0):
                break
        # print(self.reported_cells)
        self.alignment = [new_v[::-1], new_w[::-1]]

    def __global_align(self):
        """
        Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
        computed by traceback_global.

        :param: v
        :param: w
        :param: delta
        """
        # YOUR CODE HERE
        # i : length of string v and j: length of string w
        # Initialisation - score for trivial cases
        
        m = len(self.v)
        n = len(self.w)
        for i in range(0, m+1):
            for j in range (0, n+1):
                if(i==0 and j==0):
                    self.dp_table[i][j] = 0
                    self.pointers[i][j] = ORIGIN
                elif(i==0 and j>0):
                    self.dp_table[i][j] = self.dp_table[i][j-1] + self.score_func["-"][self.w[j-1]]
                    self.pointers[i][j] = LEFT
                elif(i>0 and j==0):
                    self.dp_table[i][j] = self.dp_table[i-1][j] + self.score_func[self.v[i-1]]["-"]
                    self.pointers[i][j] = UP
                else:
                    match_or_mismatch = self.dp_table[i-1][j-1] + self.score_func[self.v[i-1]][self.w[j-1]]
                    deletion = self.dp_table[i-1][j] + self.score_func[self.v[i-1]]["-"]
                    insertion = self.dp_table[i][j-1] + self.score_func["-"][self.w[j-1]]
                    self.dp_table[i][j] = max(match_or_mismatch, deletion, insertion)

                    # Update pointers
                    if self.dp_table[i][j] == match_or_mismatch:
                        self.pointers[i][j] = TOPLEFT
                    elif self.dp_table[i][j] == deletion:
                        self.pointers[i][j] = UP
                    else:
                        self.pointers[i][j] = LEFT
             
        # find the alignment by traceback on pointers table
        self.__traceback_global()

        self.global_score = self.dp_table[m][n]

    def get_alignment(self):
        start = time.time()
        self.__global_align()
        end = time.time()
        self.alignment_time = end - start
        return ''.join(self.alignment[0])+'\t'+ ''.join(self.alignment[1])
    
    def get_alignment_score(self):
        start = time.time()
        score = NeedlemanWunsch.compute_score_quad_space(self.v, self.w, self.score_func)
        end = time.time()
        self.score_time = end - start
        return score
    
    def draw_dp_table(self, ax=[], cell_size=40, score_size=16, op_file="nw_dp.png"):
        color_cell = "lightgray"
        color_path = "tomato"
        color_mid = "red"
        color_txt = "black"	

        # horizontal string
        ax.plot(0, len(v) + 1, markersize=score_size, marker="$0$", c="black")
        for ind_w in range(len(w)):
            ax.plot(ind_w + 1, len(v) + 1, markersize=score_size, marker="$" + w[ind_w] + "$", c="black")
        # vertical string
        ax.plot(-1, len(v), markersize=score_size, marker="$0$", c="black")
        for ind_v in range(len(v)):
            ax.plot(-1, ind_v, markersize=score_size, marker="$" + v[-ind_v-1] + "$", c="black")

        for ind_row in range(len(self.dp_table)):
            for ind_col in range(len(self.dp_table[0])):
                # print(f"lhs: {len(self.dp_table)-ind_row-1}, rhs: {ind_col}")
                if (len(self.dp_table)-ind_row-1, ind_col) in self.reported_cells:
                    # print("Path")
                    ax.plot(ind_col, ind_row, markersize=cell_size, marker="s", c=color_path, alpha=0.5)
                else:
                    # print("Not path")
                    ax.plot(ind_col, ind_row, markersize=cell_size, marker="s", c=color_cell, alpha=0.5)
                ax.plot(ind_col, ind_row, markersize=score_size, marker="$" + str(self.dp_table[-ind_row-1][ind_col]) + "$", c=color_txt)			
        plt.title("Alignment score dynamic programming matrix from Needleman Wunsch Algorithm")
        plt.savefig(op_file)      

if __name__ == "__main__":

    # v = "TAGATA"
    # w = "GTAGGCTTAAGGTTA"

    # v = "ATGTC"
    # w = "ATCGC"
    v = "TG"
    w = "ATCG"
    nw_obj = NeedlemanWunsch(v,w)

    alignment = nw_obj.get_alignment()
    print(f"v = {v}")
    print(f"w = {w}")
    print(f"Optimal Alignment score: {nw_obj.get_alignment_score()}")
    print(alignment)
    print(nw_obj.reported_cells)

    cell_size = 50/(1+0.04*max(len(v), len(w)))
    score_size = 30/(1+0.04*max(len(v), len(w)))
    fig, ax = plt.subplots(figsize=(5*np.sqrt(len(w)), 5*np.sqrt(len(v))))
    nw_obj.draw_dp_table(ax=ax,cell_size=cell_size, score_size=score_size)
