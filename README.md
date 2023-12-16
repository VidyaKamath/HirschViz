# HirschViz
Hirschberg Linear Space Global Alignment Visualizer


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

# Quick Start
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
- This command will generate multiple *.png files and two *.gif files. For sample output files, see results/visualization
- `recursion_tree.png`: A static image file showing the recursion tree with all the intermediate variable values.
- `recursion_tre.gif`: Animation showing the order of recursive function calls for Hirschberg Algorithm
  ![recurstion_tree.gif](https://github.com/VidyaKamath/HirschViz/blob/main/results/visualization/recursion_tree.gif)
- `dp_recur.gif`: Animation showing the recursive Dynamic programming for optimal alignment path search.
  ![dp_recur.gif](https://github.com/VidyaKamath/HirschViz/blob/main/results/visualization/dp_recur.gif)


