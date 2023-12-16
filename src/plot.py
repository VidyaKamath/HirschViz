from matplotlib import pyplot as plt
import matplotlib
import pandas as pd
import os

current_script_dir = os.path.dirname(os.path.abspath(__file__))

score_bmk_csv = os.path.join(current_script_dir, "..","results", "evaluation", "score_benchmark_0_8000.csv")
alignment_bmk_csv = os.path.join(current_script_dir, "..", "results", "evaluation", "alignment_benchmark_0_8000.csv")
fig_folder = os.path.join(current_script_dir, "..","results", "evaluation")

matplotlib.rcParams['figure.figsize'] = (40, 20)  # Increase the size of the plot
matplotlib.rcParams['font.family'] = 'sans-serif'  # Set the font to Times New Roman
matplotlib.rcParams['font.size'] = 30
score_bmk = pd.read_csv(score_bmk_csv)
alignment_bmk = pd.read_csv(alignment_bmk_csv)

# Plot space - score
xs = score_bmk["seq_length"]
ynws  = score_bmk["nw_max_mem"]
yhbs = score_bmk["hb_max_mem"]

plt.subplot(1,2,1)
plt.plot(xs,yhbs,label = 'Hirschberg', color="red")
plt.plot(xs,ynws,label = 'Needleman-Wunsch', color="blue")
plt.xlabel('Sequence Length')
plt.ylabel('Max RAM usage in Mib')
plt.title("Space analysis for computing alignment score")
plt.tight_layout(pad=2.5)
plt.legend()
plt.grid()


# Plot space - alignment
xa = alignment_bmk["seq_length"]
ynwa  = alignment_bmk["nw_max_mem"]
yhba = alignment_bmk["hb_max_mem"]

plt.subplot(1,2,2)
plt.plot(xa,yhba,label = 'Hirschberg', color="red")
plt.plot(xa,ynwa,label = 'Needleman-Wunsch', color="blue")
plt.xlabel('Sequence Length')
plt.ylabel('Max RAM usage in Mib')
plt.title("Space analysis for finding the alignment")
plt.tight_layout(pad=2.5)
plt.legend()
plt.grid()

plt.savefig(os.path.join(fig_folder, "space_analysis.png"))

plt.cla()
plt.clf()

# Plot time - score

xs = score_bmk["seq_length"]
ynws  = score_bmk["nw_time"]
yhbs = score_bmk["hb_time"]

plt.subplot(1,2,1)
plt.plot(xs,yhbs,label = 'Hirschberg', color="red")
plt.plot(xs,ynws,label = 'Needleman-Wunsch', color="blue")
plt.xlabel('sequence length')
plt.ylabel('time in seconds')
plt.title("Execution time for computing alignment score")
plt.tight_layout(pad=2.5)
plt.legend()
plt.grid()

# Plot time - alignment
xa = alignment_bmk["seq_length"]
ynwa  = alignment_bmk["nw_time"]
yhba = alignment_bmk["hb_time"]

plt.subplot(1,2,2)
plt.plot(xa,yhba,label = 'Hirschberg', color="red")
plt.plot(xa,ynwa,label = 'Needleman-Wunsch', color="blue")
plt.xlabel('sequence length')
plt.ylabel('time in seconds')
plt.title("Execution time for finding the alignment")
plt.tight_layout(pad=2.5)
plt.legend()
plt.grid()
plt.savefig(os.path.join(fig_folder, "exec_time_analysis.png"))

