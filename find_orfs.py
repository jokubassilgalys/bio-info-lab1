import os
import sys
import argparse
import time
import numpy as np
from Bio import SeqIO
from itertools import product
from collections import Counter

START_CODONS = ["ATG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]

CODONTAB = {
    'TCA': 'S',    
    'TCC': 'S',    
    'TCG': 'S',    
    'TCT': 'S',    
    'TTC': 'F',    
    'TTT': 'F',    
    'TTA': 'L',    
    'TTG': 'L',    
    'TAC': 'Y',    
    'TAT': 'Y',    
    'TAA': '*',    
    'TAG': '*',    
    'TGC': 'C',    
    'TGT': 'C',    
    'TGA': '*',    
    'TGG': 'W',    
    'CTA': 'L',    
    'CTC': 'L',    
    'CTG': 'L',    
    'CTT': 'L',    
    'CCA': 'P',    
    'CCC': 'P',    
    'CCG': 'P',    
    'CCT': 'P',    
    'CAC': 'H',    
    'CAT': 'H',    
    'CAA': 'Q',    
    'CAG': 'Q',    
    'CGA': 'R',    
    'CGC': 'R',    
    'CGG': 'R',    
    'CGT': 'R',    
    'ATA': 'I',    
    'ATC': 'I',    
    'ATT': 'I',    
    'ATG': 'M',    
    'ACA': 'T',    
    'ACC': 'T',    
    'ACG': 'T',    
    'ACT': 'T',    
    'AAC': 'N',    
    'AAT': 'N',    
    'AAA': 'K',    
    'AAG': 'K',    
    'AGC': 'S',    
    'AGT': 'S',    
    'AGA': 'R',    
    'AGG': 'R',    
    'GTA': 'V',    
    'GTC': 'V',    
    'GTG': 'V',    
    'GTT': 'V',    
    'GCA': 'A',    
    'GCC': 'A',    
    'GCG': 'A',    
    'GCT': 'A',    
    'GAC': 'D',    
    'GAT': 'D',    
    'GAA': 'E',    
    'GAG': 'E',    
    'GGA': 'G',    
    'GGC': 'G',    
    'GGG': 'G',    
    'GGT': 'G'     
}

def main(path, method):

    if not os.path.exists(path):
        print(f"Path does not exist: {path}")
        sys.exit(1)
    if os.path.isfile(path):
        print(f"Input is a file: {path}")
        process_fasta_file(path, method)
    elif os.path.isdir(path):
        print(f"Input is a directory: {path}")
        for filename in os.listdir(path):
            if filename.endswith(".fasta"):
                process_fasta_file(os.path.join(path, filename), method)
    else:
        print(f"Input is neither a file nor a directory: {path}")
        sys.exit(1)

def process_fasta_file(fasta_file, method, show_examples=5):
    print(f"\n=== Processing file: {fasta_file} using method: {method} ===")
    try:
        records = list(SeqIO.parse(fasta_file, "fasta"))
    except Exception as e:
        print(f"ERROR reading {fasta_file}: {e}")
        return

    if not records:
        print(f"No FASTA records found in {fasta_file}.")
        return


    for record in records:
        seq = str(record.seq)
        print(f"\nRecord: {record.id}")
        print(f"Length: {len(seq)} bp")

        if method == "original":
            forward_orfs = find_orfs(seq)
            rev_seq = reverse_complement(seq)
            reverse_orfs = find_orfs(rev_seq)
        elif method == "alternative":
            forward_orfs = find_orfs_alternative(seq)
            rev_seq = reverse_complement(seq)
            reverse_orfs = find_orfs_alternative(rev_seq)

        print(f"Found {len(forward_orfs)} ORFs on forward strand")
        print(f"Found {len(reverse_orfs)} ORFs on reverse complement strand")

        forward_orfs = filter_orfs(forward_orfs, 100)
        reverse_orfs = filter_orfs(reverse_orfs, 100)

        print(f"Found {len(forward_orfs)} ORFs (>=100bp) on forward strand")
        print(f"Found {len(reverse_orfs)} ORFs (>=100bp) on reverse complement strand\n")

        for orf in forward_orfs:
            orf['protein'] = translate_orf(orf['sequence'])
        for orf in reverse_orfs:
            orf['protein'] = translate_orf(orf['sequence'])

        for strand, orfs in [("Forward", forward_orfs), ("Reverse", reverse_orfs)]:
            print(f"--- {strand} strand ORFs (first {show_examples}) ---")
            if not orfs:
                print("  (none)")
            for orf in orfs[:show_examples]:
                print(f"  Frame {orf['frame']}: {orf['start']}-{orf['end']} ({orf['length']} bp)")
                print(f"  Protein (len {len(orf['protein'])} aa): {orf['protein']}")

        compare_all_sequences([orf['protein'] for orf in forward_orfs + reverse_orfs])

    print("")

def filter_orfs(orfs, min_length=100):
    return [orf for orf in orfs if orf['length'] >= min_length]
    
def reverse_complement(seq):
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]

### Translation to proteins ### 

def translate_codon(codon):
    codon = codon.upper()
    return CODONTAB.get(codon)

def translate_orf(orf_seq):
    protein = []
    seq_len = len(orf_seq)
    for i in range(0, seq_len - 2, 3):
        codon = orf_seq[i:i+3]
        aa = translate_codon(codon)
        if aa is None:
            continue
        protein.append(aa)
    return ''.join(protein)

###

### Codon and Dicodon Frequencies and Distance Matrix ###

def codon_frequencies(protein_seq):
    codons = list("ACDEFGHIKLMNPQRSTVWY*")  # 20 standard + stop
    freq = Counter(protein_seq)
    return {a: freq.get(a, 0) for a in codons}

def dicodon_frequencies(protein_seq):
    dicodons = [protein_seq[i:i+2] for i in range(len(protein_seq)-1) if len(protein_seq[i:i+2]) == 2]
    all_dicodons = [''.join(p) for p in product(list("ACDEFGHIKLMNPQRSTVWY*"), repeat=2)]
    freq = Counter(dicodons)
    return {d: freq.get(d, 0) for d in all_dicodons}

def distance_matrix(frequency_list, metric='euclidean'):
    # checking if any input frequency list is empty or inconsistent
    if not frequency_list or any(len(f) == 0 for f in frequency_list):
        raise ValueError("Empty or invalid frequency list passed to distance_matrix.")

    # converting all to arrays with consistent key order
    keys = list(frequency_list[0].keys())
    freq_arrays = [np.array([f[k] for k in keys]) for f in frequency_list]

    n = len(freq_arrays)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if metric == 'euclidean':
                matrix[i, j] = np.linalg.norm(freq_arrays[i] - freq_arrays[j])
            elif metric == 'cosine':
                a, b = freq_arrays[i], freq_arrays[j]
                matrix[i, j] = 1 - np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    return matrix

def compare_all_sequences(protein_sequences):
    codons_freqs = [codon_frequencies(p) for p in protein_sequences]
    dicodon_freqs = [dicodon_frequencies(p) for p in protein_sequences]

    # Print codons and dicodons with frequency > 0 for each sequence
    for idx, seq in enumerate(protein_sequences):
        print(f"\nSequence {idx+1} (length {len(seq)} aa):")

        codon_freq = codons_freqs[idx]
        dicodon_freq = dicodon_freqs[idx]

        codon_line = ' '.join([f"{codon}:{count}" for codon, count in codon_freq.items() if count > 0])
        print("Codons with frequency > 0:", codon_line)

        dicodon_line = ' '.join([f"{d}:{count}" for d, count in dicodon_freq.items() if count > 0])
        print("Dicodons with frequency > 0:", dicodon_line)

    codon_dist = distance_matrix(codons_freqs, 'euclidean')
    dicodon_dist = distance_matrix(dicodon_freqs, 'euclidean')

    print("\nCodon Distance Matrix:\n", codon_dist)
    print("Dicodon Distance Matrix:\n", dicodon_dist)

    return codon_dist, dicodon_dist

###

#
# --method original
# active by default
# finds codon pairs by finding the stops first and backtracking to the farthest start
#
def find_orfs(seq):
    orfs = []
    seq_len = len(seq)

    for frame in range(3):  # reading frame offset 0, 1, 2
        codons = [seq[i:i+3] for i in range(frame, seq_len-2, 3)]   # list of codons in this frame
        codon_positions = list(range(frame, seq_len-2, 3))

        for j, codon in enumerate(codons):
            if codon in STOP_CODONS:    # checking all codons for a stop
                farthest_start_index = None
                for k in range(j-1, -1, -1):    # backtracking until a stop memorising the farthest start
                    if codons[k] in STOP_CODONS:
                        break
                    if codons[k] in START_CODONS:
                        farthest_start_index = k
                if farthest_start_index is not None:    # if any start was found, append the ORF list
                    orf_seq = seq[codon_positions[farthest_start_index]:codon_positions[j]+3]
                    orfs.append({
                        "frame": frame + 1,
                        "start": codon_positions[farthest_start_index] + 1,
                        "end": codon_positions[j] + 3,
                        "length": len(orf_seq),
                        "sequence": orf_seq
                    })
    return orfs

#
# --method alternative
# finds codon pairs without backtracking by matching earliest start to the first stop
#
def find_orfs_alternative(seq):
    orfs = []
    seq = seq.upper()
    seq_len = len(seq)

    for frame in range(3):  # reading frame offset 0, 1, 2
        start_positions = [i for i in range(frame, seq_len-2, 3) if seq[i:i+3] in START_CODONS] # list of all start codons
        stop_positions = [i for i in range(frame, seq_len-2, 3) if seq[i:i+3] in STOP_CODONS]   # list of all stop codons

        combined = [(pos, 'start') for pos in start_positions] + [(pos, 'stop') for pos in stop_positions]
        combined.sort(key=lambda x: x[0])   # list of only starts and stops

        last_start = None
        for pos, codon_type in combined:
            if codon_type == 'start' and last_start is None: # if start without an earlier unmatched start
                last_start = pos
            elif codon_type == 'stop' and last_start is not None: # if found a stop with a memorised start
                orf_seq = seq[last_start:pos+3]
                orfs.append({
                    'frame': frame + 1,
                    'start': last_start + 1,
                    'end': pos + 3,
                    'length': len(orf_seq),
                    'sequence': orf_seq
                })
                last_start = None
    return orfs

#
# --compare
# optional statistics to compare find_orfs and find_orfs_alternative
#
def compare_methods(fasta_file, show_examples=5):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for record in records:
        seq = str(record.seq)
        print(f"\nRecord: {record.id}")
        
        start_time = time.time()
        orfs1 = find_orfs(seq)
        time1 = time.time() - start_time
        print(f"Method 1: {len(orfs1)} ORFs found in {time1:.6f} seconds")

        start_time = time.time()
        orfs2 = find_orfs_alternative(seq)
        time2 = time.time() - start_time
        print(f"Method 2: {len(orfs2)} ORFs found in {time2:.6f} seconds")

        print(f"Time difference: {abs(time1 - time2):.6f} seconds")

        set1 = set((orf['start'], orf['end'], orf['frame']) for orf in orfs1)
        set2 = set((orf['start'], orf['end'], orf['frame']) for orf in orfs2)
        common = set1 & set2
        only1 = set1 - set2
        only2 = set2 - set1

        print(f"Common ORFs: {len(common)}")
        print(f"Only in Method 1: {len(only1)}")
        print(f"Only in Method 2: {len(only2)}")

        if show_examples > 0:
            print("Examples from Method 1:")
            for orf in orfs1[:show_examples]:
                print(f"  Frame {orf['frame']}: {orf['start']}-{orf['end']} ({orf['length']} bp)")
            print("Examples from Method 2:")
            for orf in orfs2[:show_examples]:
                print(f"  Frame {orf['frame']}: {orf['start']}-{orf['end']} ({orf['length']} bp)")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run main or comparison mode for ORF finding.")
    parser.add_argument("fasta_file", help="Path to a FASTA file.")
    parser.add_argument("--compare", action="store_true", help="Run comparison between both ORF finding methods.")
    parser.add_argument("--method", choices=["original", "alternative"], default="original", help="Choose which ORF finding method to use in main processing.")
    args = parser.parse_args()

    if args.compare:
        compare_methods(args.fasta_file)
    else:
        main(args.fasta_file, method=args.method)