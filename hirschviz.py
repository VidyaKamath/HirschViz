import argparse
from src.hirschberg_viz import hirsch_viz_main

def main():
    parser = argparse.ArgumentParser(description='')

    # Required string arguments
    parser.add_argument('-o', type=str, required=True,help="Output directory to save visualization files")
    parser.add_argument('-v', type=str, required=True, help='Sequence 1')
    parser.add_argument('-w', type=str, required=True, help='Sequence 2')

    args = parser.parse_args()

    # Accessing the provided arguments
    v = args.v
    w = args.w
    op_dir = args.o  
    hirsch_viz_main(v,w,op_dir)
    
if __name__ == "__main__":
    main()
