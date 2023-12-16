from src.score_functions import get_blosum62, get_dna_score_func
import numpy as np
import time
class Hirschberg(object):

    def __init__(self, v, w):
        self.v = v
        self.w = w
        self.score_func = get_dna_score_func()
        
        # Output
        self.alignment = []
        self.global_score = None
        self.reported_cells = []

        self.alignment_time = None
        self.score_time = None
    
    @staticmethod
    def compute_score_lin_space(v, w, delta):
        """
        _summary_

        Computes the prefix and the suffix score using only two columns and len(v) rows
        """
        m = len(v)
        n = len(w)
        M = [[0 for j in range(2)] for i in range(m+1)]
        col = 0
        for j in range(n + 1):
            for i in range(m + 1):
                if i==0 and j==0:
                    M[i][col] = 0
                elif i==0 and j>0:
                    M[i][col] = M[i][1-col] + delta["-"][w[j-1]]
                elif i>0 and j==0:
                    M[i][col] = M[i-1][col] + delta[v[i-1]]["-"]
                else:
                    match_or_mismatch = M[i-1][1-col] + delta[v[i-1]][w[j-1]]
                    deletion = M[i-1][col] + delta[v[i-1]]["-"]
                    insertion = M[i][1-col] + delta["-"][w[j-1]]
                    M[i][col] = max(match_or_mismatch,deletion,insertion)
            col = 1-col

        last_column = [row[1-col] for row in M]
        score = M[m][1-col]
        return score, last_column
        
    # To get the desired tree from HW1 - new
    @staticmethod
    def __hirschberg(v,w,i,j,ip,jp,res,score_func):
        # base case
        if jp-j==1:
            res.append((i, j))
            res.append((ip, jp))
            return 
        if jp - j < 1:
            return
        
        # middle column to split
        mid = int((j+jp)/2)
        _, prefix = Hirschberg.compute_score_lin_space(v[i:ip:],w[j:mid:], score_func)
        _, suffix = Hirschberg.compute_score_lin_space(v[i:ip:][::-1],w[mid:jp:][::-1], score_func)
        
        maxidx = -1
        maxweight = float("-inf")

        for idx in range(len(prefix)):
            weight = prefix[idx] + suffix[-idx-1]
            if maxidx == -1 or maxweight <= weight:
                maxidx = idx
                maxweight = weight
        
        # Report the cell on the traceback path
        res.append((maxidx+i, mid))
        
        # print("recursive call ",cc)
        # print(f"i*: {maxidx+i}, Prefix: {prefix[maxidx]}, Suffix: {suffix[-maxidx-1]}, report: ({maxidx+i}, {mid})")
        # print(f"Left i:{i},j:{j}, ip:{maxidx+i},jp:{mid}")
        # print(f"Left {cc} ({i},{j},{maxidx+i},{mid})")
        Hirschberg.__hirschberg(v=v, w=w, i=i, j=j, ip=maxidx+i, jp=mid, res=res, score_func=score_func)
        # print(f"Right i:{maxidx+i},j:{mid}, ip:{ip},jp:{jp}")
        # print(f"Right {cc} ({maxidx+i},{mid},{ip},{jp})")
        Hirschberg.__hirschberg(v=v, w=w, i=maxidx+i, j=mid, ip=ip, jp=jp, res=res, score_func=score_func)

    @ staticmethod
    def reconstruct_path(s1, s2, cells):
        """
        _summary_
        Reference: This function was taken from a GitHub repo. Todo: Update this to use column-wise instead of row-wise
        """
        len1 = len(s1)
        len2 = len(s2)
        dp = [[0 for j in range(len1 + 1)] for i in range(2)]
        curridx = 0
        for i in range(len2 + 1):
            for j in range(len1 + 1):
                if i==0 and j==0:
                    dp[curridx][j]=0
                elif i>0 and j==0:
                    dp[curridx][j]=dp[1-curridx][j]-1
                elif i==0 and j>0:
                    dp[curridx][j]=dp[curridx][j-1]-1
                else:
                    score = 1 if s1[j-1] == s2[i-1] else -1
                    dp[curridx][j] = max(dp[1-curridx][j]-1,max(dp[curridx][j-1]-1,dp[1-curridx][j-1]+score))
            curridx = 1-curridx
        
            left = curridx
            right = 1 - curridx
            
            # fill in missing cells
            if i > 0 and cells[i][0] - cells[i - 1][0] > 1:
                start_cell = cells[i - 1]
                start_row, start_col = start_cell

                end_cell = cells[i]
                end_row, end_col = end_cell
                
                stack = []
                stack.append(start_cell)
                bt = {}
                
                # shortest path from start cell to end cell given constraint
                while True:
                    current_cell = stack.pop()
                    cell_row, cell_col = current_cell

                    if current_cell == end_cell:
                        break

                    # go right
                    if cell_col < end_col and dp[left][cell_row] - 1 == dp[right][cell_row]:
                        next_cell = (cell_row, cell_col + 1)
                        stack.append(next_cell)
                        bt[next_cell] = current_cell
                        
                    # go down
                    column = right if cell_col == end_col else left
                    if cell_row < end_row and dp[column][cell_row] - 1 == dp[column][cell_row + 1]:
                        next_cell = (cell_row + 1, cell_col)
                        stack.append(next_cell)
                        bt[next_cell] = current_cell
                        
                    # go diagonal - match
                    if cell_row < end_row and cell_col < end_col and s1[cell_row] == s2[cell_col] and dp[left][cell_row] + 1 == dp[right][cell_row + 1]:
                        new_cell = (cell_row + 1, cell_col + 1)
                        stack.append(new_cell)
                        bt[new_cell] = current_cell
                        
                    # go diagonal - mismatch
                    if cell_row < end_row and cell_col < end_col and s1[cell_row] != s2[cell_col] and dp[left][cell_row] - 1 == dp[right][cell_row + 1]:
                        new_cell = (cell_row + 1, cell_col + 1)
                        stack.append(new_cell)
                        bt[new_cell] = current_cell
                        
                # backtrace
                current_cell = bt[end_cell]
                while current_cell != start_cell:
                    cells.append(current_cell)
                    current_cell = bt[current_cell]
        return sorted(cells)
    
    def reconstruct_alignment(self):
    
        seq1 = ""
        seq2 = ""

        for s in self.reported_cells:
            if (s == (0,0)):
                temp = s
                continue
            else:
                if(s[0]-temp[0] == 0):
                    seq1 += "-"
                    seq2 += self.w[s[1]-1]
                elif(s[0]-temp[0] == 1):
                    seq1 += self.v[s[0]-1]
                    seq2 += self.w[s[1]-1]
                else:
                    seq1 += self.v[temp[0]]
                    seq2 += self.w[temp[1]]
                    for j in range(s[0]-temp[0]-1):
                        seq2 += "-"
                        seq1 += self.v[temp[0]+1+j]
                temp = s
        self.alignment = [seq1, seq2]

    def get_alignment(self):
        start = time.time()
        res = []
        m = len(self.v)
        n = len(self.w)
        self.__hirschberg(v=self.v, w=self.w, i=0, j=0, ip=m, jp=n, res=res, score_func=self.score_func)
        res = list(set(res))
        self.reported_cells = sorted(res,key=lambda x: (x[0],x[1]))

        # reconstruct path - ONLY FOR VISUALISATION
        # self.reported_cells = Hirschberg.reconstruct_path(self.v, self.w, path)

        # construct the alignment
        self.reconstruct_alignment()
        end = time.time()
        self.alignment_time = end - start
        return ''.join(self.alignment[0])+'\t'+ ''.join(self.alignment[1])
    
    def get_alignment_score(self):
        start = time.time()
        score, _ = Hirschberg.compute_score_lin_space(self.v, self.w, self.score_func)
        end = time.time()
        self.score_time = end-start
        return score

    
if __name__ == "__main__":
    
    # test 1
    v = "TAGATA"
    w = "GTAGGCTTAAGGTTA"

    # v = "ATGTC"
    # w = "ATCGC"

    hb_obj = Hirschberg(v,w)

    alignment = hb_obj.get_alignment()
    print(f"v = {v}")
    print(f"w = {w}")
    print(f"Optimal Alignment score: {hb_obj.get_alignment_score()}")
    print(alignment)
    print(hb_obj.reported_cells)
