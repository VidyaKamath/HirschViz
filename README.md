# HirschViz
Hirschberg Linear Space Global Alignment Visualizer

There are many interactive tools available for visualizing the Needleman-Wunsch algorithm like [Alignment Visualizer](https://valiec.github.io/AlignmentVisualizer/index.html}) and [Interactive NW](http://experiments.mostafa.io/public/needleman-wunsch/). However, in the case of Hirschberg's algorithm, a quick Google search leads only to some online content showcasing Hirschberg's visualization through YouTube videos [Hirschebrg Video](https://www.youtube.com/watch?v=cPQeJt-2Y1Q) and static images. Most visualizations only show the final result without detailing the step-by-step process, especially the recursive steps that are crucial for a deeper understanding. Also, a quick online search did not yield any tool offering a detailed, step-by-step view of how Hirschberg's algorithm fills its tables recursively. This gap in available resources motivates the need to create this `HirschViz` tool that focuses on showing these recursive steps in a clear and accessible way, aiming to help students better grasp the internals of the algorithm.

# Dependencies
- Python installation version > 3.9 (This project uses version 3.9.6)
- A virtual environment can be created using python `virtualenv`

# Setup
## Clone the repository
```
git clone https://github.com/VidyaKamath/HirschViz.git
```
## Create a virtual environment
```
pip install virtualenv
python3 -m venv <env_path>
```
## Activate the virtual environment
```
source <env_path>/bin/activate
```
## Install the necessary packages
```
pip install -r requirements.txt
```
# Available Commands
## To generate step-by-step visualization for Hirschberg Algorithm use the below command. 
- See [Quick Start](#Quick_Start) for an example command.
```
(hirsch_viz) vkpailodi@Vidyas-MacBook-Air HirschViz % python3 hirschviz.py -h
usage: hirschviz.py [-h] -o O -v V -w W

optional arguments:
  -h, --help  show this help message and exit
  -o O        Output directory to save visualization files
  -v V        Sequence 1
  -w W        Sequence 2
```
## To run benchmark use the below command 
- See [Quick Start](#running-benchmark) for an example command.
```
(hirsch_viz) vkpailodi@Vidyas-MacBook-Air HirschViz % python3 run_benchmark.py -h
usage: run_benchmark.py [-h] [-o O] --start START --stop STOP --step STEP [--score-only] [--both] 

Argument Parser Example

optional arguments:
  -h, --help     show this help message and exit
  -o O           csv file to save results
  --start START  start sequence length
  --stop STOP    stop sequence length
  --step STEP    step size
  --score-only   Run benchmark to compute score / otherwise to find the alignment
  --both         Run benchmark to compute score and to find the alignment
```
# Quick_Start
- Make sure the current directory is the root folder of the cloned repository directory
- Make sure the Setup is successful and the environment is activated
## Run Hirschberg Visualization
```
python3 hirschviz.py -o ~/hirsch_vizop/ -v TG -w ATCG
```
The console output should be this
  ```
  (hirsch_viz) vkpailodi@Vidyas-MacBook-Air HirschViz % python3 hirschviz.py -o ~/hirsch_vizop/ -v TG -w ATCG
  Starting to make animation
  File /Users/vkpailodi/hirsch_vizop/recursion_tree.png successfully written
  Writing frames....
  Writing gif...
  Saved gif /Users/vkpailodi/hirsch_vizop/recursion_tree.gif successfully
  ```
- This command will generate multiple *.png files and two *.gif files. For sample output files, see [results](https://github.com/VidyaKamath/HirschViz/blob/main/results/visualization/)
- `recursion_tree.png`: A static image file showing the recursion tree with all the intermediate variable values.
  ![recurstion_tree.png](https://github.com/VidyaKamath/HirschViz/blob/main/results/visualization/recursion_tree.png)
- `recursion_tre.gif`: Animation showing the order of recursive function calls for Hirschberg Algorithm
  ![recurstion_tree.gif](https://github.com/VidyaKamath/HirschViz/blob/main/results/visualization/recursion_tree.gif)
- `dp_recur.gif`: Animation showing the recursive Dynamic programming for optimal alignment path search.
  ![dp_recur.gif](https://github.com/VidyaKamath/HirschViz/blob/main/results/visualization/dp_recur.gif)

## Running Benchmark
### Computing optimal alignment score: linear space v/s quadratic space
```
python3 run_benchmark.py --start 0 --stop 8500 --step 500 -o score_benchmark.csv --score-only
```
Output file: [score_benchmark.csv](https://github.com/VidyaKamath/HirschViz/blob/main/results/evaluation/score_benchmark_0_8000.csv)

### Finding optimal alignment: linear space v/s quadratic space
```
python3 run_benchmark.py --start 0 --stop 8500 --step 500 -o alignment_benchmark.csv
```
Output file: [alignment_benchmark.csv](https://github.com/VidyaKamath/HirschViz/blob/main/results/evaluation/alignment_benchmark_0_8000.csv)

### Computing optimal score and  alignment: linear space v/s quadratic space
```
python3 run_benchmark.py --start 0 --stop 10 --step 1 -o score_alignment_benchmark.csv --both
```
Output file: [score_alignment_benchmark.csv](https://github.com/VidyaKamath/HirschViz/blob/main/results/evaluation/correctness_check_alignment.csv)


 


