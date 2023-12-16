from src.score_functions import get_dna_score_func
import numpy as np
import os
from visualiser.visualiser import Visualiser as vs
import matplotlib.pyplot as plt
import imageio.v2 as imageio


class HirschbergViz(object):
    level_dict = {}
    res = []

    def __init__(self, v, w, op_dir):
        self.v = v
        self.w = w
        self.score_func = get_dna_score_func()
        self.cell_size = 50/(1+0.04*max(len(v), len(w)))
        self.score_size = 30/(1+0.04*max(len(v), len(w)))

        # Output
        self.alignment = []
        self.global_score = None
        self.reported_cells = []

        self.alignment_time = None
        self.score_time = None

        self.op_dir = op_dir
    
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
    
    @staticmethod
    def global_align(v, w, delta):
        """
        Returns the score of the maximum scoring alignment of the strings v and w, as well as the actual alignment as
        computed by traceback_global.
        :param: v
        :param: w
        :param: delta
        """
        M = [[0 for j in range(len(w)+1)] for i in range(len(v)+1)]
        score = None
        for i in range(0, len(v)+1):
            for j in range (0, len(w)+1):
                if(i==0 and j==0):
                    M[i][j] = 0
                elif(i==0 and j>0):
                    M[i][j] = M[i][j-1] + delta["-"][w[j-1]]
                elif(i>0 and j==0):
                    M[i][j] = M[i-1][j] + delta[v[i-1]]["-"]
                
        m = len(v)
        n = len(w)
        for i in range(1, m+1):
            for j in range(1, n+1):
                match_or_mismatch = M[i-1][j-1] + delta[v[i-1]][w[j-1]]
                deletion = M[i-1][j] + delta[v[i-1]]["-"]
                insertion = M[i][j-1] + delta["-"][w[j-1]]
                M[i][j] = max(match_or_mismatch, deletion, insertion)
        score = M[m][n]
        return score, M
    
    @vs(ignore_args=["res","V","W","score_func"], node_properties_kwargs={"shape":"record", "color":"#f57542", "style":"filled", "fillcolor":"grey", "imagescale":True})
    def __hirschberg_viz_recur_tree(v,w,i,j,ip,jp,res,score_func,level):
        if jp-j==1:
            res.append((i, j))
            res.append((ip, jp))
            return 
        if jp - j < 1:
            return
        
        mid = int((j+jp)/2)
        # print(f"Prefix score: {s1[i:ip:]} {s2[j:mid:]}")
        p_score, prefix = HirschbergViz.compute_score_lin_space(v[i:ip:],w[j:mid:], score_func)
        # print(f"Suffix score: {s1[i:ip:][::-1]} {s2[mid:jp:][::-1]}")
        s_score, suffix = HirschbergViz.compute_score_lin_space(v[i:ip:][::-1],w[mid:jp:][::-1], score_func)
        # print("Global Alignemnt Prefix score")
        gp_score, dp_table_p = HirschbergViz.global_align(v[i:ip:],w[j:mid:], score_func)
        # print("Global Alignemnt Suffix score")
        gs_score, dp_table_s = HirschbergViz.global_align(v[i:ip:][::-1],w[mid:jp:][::-1], score_func)
        assert p_score == gp_score
        assert s_score == gs_score
        maxidx = -1
        maxweight = -1
        for idx in range(len(prefix)):
            weight = prefix[idx] + suffix[-idx-1]
            if maxidx == -1 or maxweight <= weight:
                maxidx = idx
                maxweight = weight
        res.append((maxidx+i, mid))
        
        # print("recursive call ",cc)
        # print(f"i*: {maxidx+i}, Prefix: {prefix[maxidx]}, Suffix: {suffix[-maxidx-1]}, report: ({maxidx+i}, {mid})")
        # print(f"Left i:{i},j:{j}, ip:{maxidx+i},jp:{mid}")
        # print(f"Left {cc} ({i},{j},{maxidx+i},{mid})")
        HirschbergViz.__hirschberg_viz_recur_tree(v=v,w=w,i=i,j=j,ip=maxidx+i,jp=mid,res=res, level=level+1, score_func=score_func)
        
        # print(f"Right i:{maxidx+i},j:{mid}, ip:{ip},jp:{jp}")
        # print(f"Right {cc} ({maxidx+i},{mid},{ip},{jp})")
        HirschbergViz.__hirschberg_viz_recur_tree(v=v,w=w,i=maxidx+i,j=mid,ip=ip,jp=jp,res=res, level=level+1, score_func=score_func)
        

         # Data for visualization
        level_data = (dp_table_p, dp_table_s, i, j, ip, jp, maxidx+i, mid)
        HirschbergViz.level_dict[level].append(level_data)
        
        rev_suffix = suffix[::-1]
        op_str = f""" 
                prefix(i) = {prefix}
                suffix(i) = {rev_suffix}\n
                (i_star, j_mid) = ({maxidx}, {mid}) \n
                prefix(i_star) = {prefix[maxidx]}
                suffix(i_star) = {rev_suffix[maxidx]}\n
                Report (i_star + i, mid): ({maxidx+i}, {mid}) \n
                Left child Input (prefix): \t
                v = {v[i:ip:]}, w = {w[j:mid:]}\n
                Right child Input (suffix): \t
                v = {v[i:ip:][::-1]}, w = {w[mid:jp:][::-1]}\n
                """
        return op_str

    @vs(ignore_args=["score_func", "op_prfx"],node_properties_kwargs={"shape":"record", "color":"#f57542", "style":"filled", "fillcolor":"grey", "imagescale":True})
    def hirschberg_viz_recur_tree_main(v,w, score_func):
        res = []
        len1 = len(v)
        len2 = len(w)
        num_levels = int(np.ceil(np.log2(len(w))))+1
        
        for ind_level in range(num_levels):
            HirschbergViz.level_dict[ind_level] = []
        # i_star, mid, prefix, suffix = get_initial_i_star_mid_and_scores(v, w,0,0,len1,len2)
        HirschbergViz.__hirschberg_viz_recur_tree(v=v,w=w,i=0,j=0,ip=len1,jp=len2, res=res, score_func=score_func,level=0)
        res = list(set(res))
        # res.sort(key=lambda x: (x[0],x[1]))
        path = sorted(res,key=lambda x: (x[0],x[1]))
        return f"Alignment:{path}", path

    @ staticmethod
    def reconstruct_path(s1, s2, cells):
        """
        _summary_

        Reference: This function was taken from a guthub repo. Update this to use column  wise
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
    

    def draw_hirschberg_dp_tbl_recur(self, level_ind=0, ax=[]):  
        v = self.v
        w = self.w

        color_cell_p = "mediumpurple"
        color_cell_s = "darkseagreen"
        color_txt = "black"    
        color_txt_p = "indigo"
        color_txt_s = "darkgreen"    
        level_data_list = HirschbergViz.level_dict[level_ind]   


        blk_sizes = [level_data_list[x][5] for x in range(len(level_data_list))]
        w_ext = "0"
        for ind_w in range(len(w)):
            w_ext = w_ext+w[ind_w]
            if ind_w+1 in blk_sizes:
                if ind_w<len(w)-1:
                    w_ext = w_ext+"00"
                else:
                    w_ext = w_ext+"0"
                
        blk_sizes = [level_data_list[x][4] for x in range(len(level_data_list))]
        v_ext = "0"
        for ind_v in range(len(v)):
            v_ext = v_ext+v[ind_v]
            if ind_v+1 in blk_sizes:
                v_ext = v_ext+"0"
    
        n_rows = len(v_ext)
        for ind_dp in range(len(level_data_list)):
            # horizontal string                        
            for ind_w in range(len(w_ext)):
                ax.plot(ind_w, len(v_ext), markersize=self.score_size, marker="$" + w_ext[ind_w] + "$", c="black")
                    
            # vertical string
            for ind_v in range(len(v_ext)-1):
                ax.plot(-1, ind_v+1, markersize=self.score_size, marker="$" + v_ext[-ind_v-2] + "$", c=color_txt_p)
            
            for ind_v in range(len(v_ext)-1):
                ax.plot(len(w_ext), ind_v+1, markersize=self.score_size, marker="$" + v_ext[-ind_v-1] + "$", c=color_txt_s)
            
            level_data = level_data_list[ind_dp]
            dp_table_p = level_data[0]
            dp_table_s = level_data[1]
            i, j, ip, jp, istar, jmid = level_data[2], level_data[3], level_data[4], level_data[5], level_data[7], level_data[7]
            
            # plot prefix table
            for ind_row in range(len(dp_table_p)):
                for ind_col in range(len(dp_table_p[0])):
                    plt_x = (2*ind_dp)+j+ind_col # col
                    plt_y = n_rows-1-(i+ind_row)-(ind_dp) # row
                    ax.plot(plt_x, plt_y, markersize=self.cell_size, marker="s", c=color_cell_p, alpha=0.5)
                    ax.plot(plt_x, plt_y, markersize=self.score_size, marker="$" + str(dp_table_p[ind_row][ind_col]) + "$", c=color_txt)

            # plot suffix table
            for ind_row in range(len(dp_table_s)):
                for ind_col in range(len(dp_table_s[0])):
                    plt_x = (2*ind_dp)+jmid+ind_col+1 # col
                    plt_y = n_rows-1-(i+ind_row)-(ind_dp) # row
                    ax.plot(plt_x, plt_y, markersize=self.cell_size, marker="s", c=color_cell_s, alpha=0.5)
                    ax.plot(plt_x, plt_y, markersize=self.score_size, marker="$" + str(dp_table_s[-ind_row-1][-ind_col-1]) + "$", c=color_txt)
        plt.savefig(os.path.join(self.op_dir, f"dp_recur_level_{level_ind}.png"))
    
    
    def draw_hirschberg_recur_leaf(self, reported_cells, dp_table, ax=[]):
        color_cell = "lightgray"
        color_path = "tomato"
        color_txt = "black"	

        # horizontal string
        ax.plot(0, len(self.v) + 1, markersize=self.score_size, marker="$0$", c="black")
        for ind_w in range(len(self.w)):
            ax.plot(ind_w + 1, len(self.v) + 1, markersize=self.score_size, marker="$" + self.w[ind_w] + "$", c="black")
        
        # vertical string
        ax.plot(-1, len(self.v), markersize=self.score_size, marker="$0$", c="black")
        for ind_v in range(len(self.v)):
            ax.plot(-1, ind_v, markersize=self.score_size, marker="$" + self.v[-ind_v-1] + "$", c="black")

        for ind_row in range(len(dp_table)):
            for ind_col in range(len(dp_table[0])):
                if (len(dp_table)-ind_row-1, ind_col) in reported_cells:
                    ax.plot(ind_col, ind_row, markersize=self.cell_size, marker="s", c=color_path, alpha=0.5)
                else:
                    ax.plot(ind_col, ind_row, markersize=self.cell_size, marker="s", c=color_cell, alpha=0.5)
                ax.plot(ind_col, ind_row, markersize=self.score_size, marker="$" + str(dp_table[-ind_row-1][ind_col]) + "$", c=color_txt)			
        
        if len(reported_cells) == 0:
            filename = "dp_recur_leaf_0.png"
        else:
            filename = "dp_recur_leaf.png"
        plt.savefig(os.path.join(self.op_dir, filename))


    def draw_dp_table(self,path):
        v = self.v
        w = self.w
        path = self.reconstruct_path(v,w, path)
        num_levels = int(np.ceil(np.log2(len(w))))
        image_files = []
        fig, ax = plt.subplots(figsize=(5*np.sqrt(len(w)), 5*np.sqrt(len(v))))
        dp_table = [["\ " for j in range(len(w)+1)] for i in range(len(v)+1)]
        self.draw_hirschberg_recur_leaf([], dp_table, ax=ax)
        image_files.append(os.path.join(self.op_dir, "dp_recur_leaf_0.png"))
        for ind in range(num_levels):	
            plt.cla()
            self.draw_hirschberg_dp_tbl_recur(level_ind = ind, ax=ax)
            file_name = f"dp_recur_level_{ind}.png"
            image_files.append(os.path.join(self.op_dir, file_name))
        plt.cla()
        score, dp_table = HirschbergViz.global_align(v, w, self.score_func)
        self.draw_hirschberg_recur_leaf(path, dp_table, ax=ax)
        image_files.append(os.path.join(self.op_dir, "dp_recur_leaf.png"))	
        return image_files


def hirsch_viz_main(v,w,op_dir):
    os.makedirs(op_dir, exist_ok=True)
    hb_obj = HirschbergViz(v,w,op_dir)

    str_path, path = HirschbergViz.hirschberg_viz_recur_tree_main(v=hb_obj.v, w=hb_obj.w, score_func=hb_obj.score_func)
    vs.make_animation(os.path.join(op_dir, "recursion_tree.gif"), delay=1000)

    # DP table recursion visualization
    images = hb_obj.draw_dp_table(path)

    # Create Gif
    gif_name = os.path.join(op_dir, "dp_recursion.gif")
    images = list(map(lambda filename: imageio.imread(filename), images))
    imageio.mimsave(gif_name, images, duration = 2000) # modify the frame duration as needed




