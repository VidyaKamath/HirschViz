def get_blosum62():
    """
    Reads in the ncbi's BLOSUM62.txt file and loads the scoring matrix
    into a dictionary.
    :return: a dictionary of dictionaries which will hold the cost of various amino acid
    substitutions as defined in BLOSUM62.
    """
    blosum = 'blosum/BLOSUM62.txt'
    delta = {}
    with open(blosum, 'r') as f:
        lines = f.readlines()[6:]
        keys = lines[0].split()
        keys[-1] = '-'
        for i, line in enumerate(lines[1:]):
            delta[keys[i]] = {k : int(v) for (k,v) in zip(keys, line.split()[1:])}
    return delta

def get_dna_score_func(match=1, other=-1):
	keys = ['A', 'C', 'T', 'G', '-']
	delta = {}
	for i in range(len(keys)):
		delta[keys[i]] = {k : v for (k,v) in zip(keys, [match if keys[i] == keys[j]  else other for j in range(len(keys))])}
	return delta