from memory_profiler import memory_usage
import random
import argparse


from src.needleman_wunsch import NeedlemanWunsch
from src.hirschberg import Hirschberg

def generate_random_DNA_sequence(length):
    bases = 'ATGC'
    random_sequence_1 = ''.join(random.choice(bases) for _ in range(length))
    random_sequence_2 = ''.join(random.choice(bases) for _ in range(length))
    return random_sequence_1, random_sequence_2


def run_score_benchmark(start, stop, step, op_csv):
    if op_csv is None:
        "score_benchmark.csv"
    headers = ["seq_length","nw_score","hb_score","nw_max_mem","hb_max_mem","nw_time","hb_time"]
    with open (op_csv, "w") as fp:
        fp.write(",".join(headers) + "\n")

        for length in range (start, stop, step):
            v, w = generate_random_DNA_sequence(length)

            # Hirschberg
            hb_obj = Hirschberg(v,w)
            # Needleman Wunsch 
            nw_obj = NeedlemanWunsch(v,w)

            mem_hb, hb_score = memory_usage(proc=(hb_obj.get_alignment_score), max_usage=True, retval=True)
            mem_nw, nw_score = memory_usage(proc=(nw_obj.get_alignment_score), max_usage=True, retval=True)

            # csv_data
            data = f"{length}, {nw_score}, {hb_score}, {mem_nw}, {mem_hb}, {nw_obj.score_time}, {hb_obj.score_time}"
            fp.write(data)
            fp.write("\n")
            print(f"Done with seq length = {length}")
    
    
def run_alignment_benchmark(start, stop, step, op_csv):
    if op_csv is None:
        "score_benchmark.csv"
    headers = ["seq_length","nw_max_mem", "hb_max_mem","nw_time", "hb_time", "nw_alignment", "hb_alignment"]

    fp = open(op_csv, "w+")
    fp.write(",".join(headers) + "\n")
    
    for length in range (start, stop, step):
        v, w = generate_random_DNA_sequence(length)

        # Hirschberg
        hb_obj = Hirschberg(v,w)
        # Needleman Wunsch 
        nw_obj = NeedlemanWunsch(v,w)

        mem_hb, alignment_hb = memory_usage(proc=(hb_obj.get_alignment), max_usage=True, retval=True)
        mem_nw, alignment_nw = memory_usage(proc=(nw_obj.get_alignment), max_usage=True, retval=True)
        data = f"{length},{mem_nw},{mem_hb},{nw_obj.alignment_time},{hb_obj.alignment_time}, {alignment_nw}, {alignment_hb}"
        fp.write(data)
        fp.write("\n")
        print(f"Done with seq length = {length}")
    fp.close()


def run_alignment_score_benchmark(start, stop, step, op_csv):
    if op_csv is None:
        "score_benchmark.csv"
    headers = ["seq_length","v","w", "nw_score", "hb_score", "nw_max_mem_a", "hb_max_mem_a","nw_time_a", "hb_time_a", "nw_max_mem_s", "hb_max_mem_s","nw_time_s", "hb_time_s", "nw_alignment", "hb_alignment"]
    fp = open(op_csv, "w+")
    fp.write(",".join(headers) + "\n")
    
    for length in range (start, stop, step):
        v, w = generate_random_DNA_sequence(length)

        # Hirschberg
        hb_obj = Hirschberg(v,w)
        # Needleman Wunsch 
        nw_obj = NeedlemanWunsch(v,w)

        mem_hbs, hb_score = memory_usage(proc=(hb_obj.get_alignment_score), max_usage=True, retval=True)
        mem_nws, nw_score = memory_usage(proc=(nw_obj.get_alignment_score), max_usage=True, retval=True)

        mem_hba, alignment_hb = memory_usage(proc=(hb_obj.get_alignment), max_usage=True, retval=True)
        mem_nwa, alignment_nw = memory_usage(proc=(nw_obj.get_alignment), max_usage=True, retval=True)
        data = f"{length},{v},{w},{nw_score}, {hb_score}, {mem_nwa},{mem_hba},{nw_obj.alignment_time},{hb_obj.alignment_time}, {mem_nws},{mem_hbs},{nw_obj.score_time}, {hb_obj.score_time}, {alignment_nw}, {alignment_hb}"
        fp.write(data)
        fp.write("\n")
        print(f"Done with seq length = {length}")
    fp.close()


def main():
    parser = argparse.ArgumentParser(description='Argument Parser Example')

    # Required string arguments
    parser.add_argument('-o', type=str, help="csv file to save results")

    # Required integer arguments
    parser.add_argument('--start', type=int, required=True, help='start Sequence length')
    parser.add_argument('--stop', type=int, required=True, help='stop sequence length')
    parser.add_argument('--step', type=int, required=True, help='step size')

    # Optional boolean argument
    parser.add_argument('--score-only', action='store_true', help='Run benchmark to compute score / otherwise to find the alignement')
    parser.add_argument('--both', action='store_true', help='Run benchmark to compute score and to find the alignement')

    args = parser.parse_args()

    # Accessing the provided arguments
    start = args.start
    stop = args.stop
    step = args.step
    score_only = args.score_only
    both = args.both
    op_csv = args.o

    # Use the arguments in your script as needed
    print(f"start: {start}, stop: {stop}, step: {step}, score_only: {score_only}")
    if score_only:
        run_score_benchmark(start,stop,step,op_csv)
    if both:
        run_alignment_score_benchmark(start,stop,step,op_csv)
    else:
        run_alignment_benchmark(start,stop,step,op_csv)

    print("Finished running benchmark")
    
if __name__ == "__main__":
    main()
