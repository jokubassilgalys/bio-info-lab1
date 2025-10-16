import os
import sys
import argparse
import time
from Bio import SeqIO

START_CODONS = ["ATG"]
STOP_CODONS = ["TAA", "TAG", "TGA"]

def main(path, method="original"):

    if not os.path.exists(path):
        print(f"Path does not exist: {path}")
        sys.exit(1)
    if os.path.isfile(path):
        print(f"Input is a file: {path}")
        process_fasta_file(path)
    elif os.path.isdir(path):
        print(f"Input is a directory: {path}")
        for filename in os.listdir(path):
            if filename.endswith(".fasta"):
                process_fasta_file(os.path.join(path, filename))
    else:
        print(f"Input is neither a file nor a directory: {path}")
        sys.exit(1)

def process_fasta_file(fasta_file, show_examples=5, method="original"):
    print(f"\n=== Processing file: {fasta_file} ===")
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

        for strand, orfs in [("Forward", forward_orfs), ("Reverse", reverse_orfs)]:
            print(f"--- {strand} strand ORFs (first {show_examples}) ---")
            if not orfs:
                print("  (none)")
            for orf in orfs[:show_examples]:
                print(f"  Frame {orf['frame']}: {orf['start']}-{orf['end']} ({orf['length']} bp)")
    print("")
     

def reverse_complement(seq):
    complement = str.maketrans("ATGC", "TACG")
    return seq.translate(complement)[::-1]

def find_orfs(seq):
    orfs = []
    seq_len = len(seq)

    # for frame in range(3):
    #     i = frame
    #     while i < seq_len - 2:
    #         codon = seq[i:i+3]
    #         if codon in START_CODONS:
    #             for j in range(i + 3, seq_len - 2, 3):
    #                 stop_codon = seq[j:j+3]
    #                 if stop_codon in STOP_CODONS:
    #                     orf_seq = seq[i:j+3]
    #                     orfs.append({
    #                         "frame": frame + 1,
    #                         "start": i + 1,
    #                         "end": j + 3,
    #                         "length": len(orf_seq),
    #                         "sequence": orf_seq
    #                     })
    #                     break
    #         i += 3
    for frame in range(3):
        codons = [seq[i:i+3] for i in range(frame, seq_len-2, 3)]
        codon_positions = list(range(frame, seq_len-2, 3))

        for j, codon in enumerate(codons):
            if codon in STOP_CODONS:
                farthest_start_index = None
                for k in range(j-1, -1, -1):
                    if codons[k] in STOP_CODONS:
                        break
                    if codons[k] in START_CODONS:
                        farthest_start_index = k
                if farthest_start_index is not None:
                    orf_seq = seq[codon_positions[farthest_start_index]:codon_positions[j]+3]
                    orfs.append({
                        "frame": frame + 1,
                        "start": codon_positions[farthest_start_index] + 1,
                        "end": codon_positions[j] + 3,
                        "length": len(orf_seq),
                        "sequence": orf_seq
                    })
    return orfs


def find_orfs_alternative(seq):
    orfs = []
    seq = seq.upper()
    seq_len = len(seq)

    for frame in range(3):
        start_positions = [i for i in range(frame, seq_len-2, 3) if seq[i:i+3] in START_CODONS]
        stop_positions = [i for i in range(frame, seq_len-2, 3) if seq[i:i+3] in STOP_CODONS]

        combined = [(pos, 'start') for pos in start_positions] + [(pos, 'stop') for pos in stop_positions]
        combined.sort(key=lambda x: x[0])

        last_start = None
        for pos, codon_type in combined:
            if codon_type == 'start':
                last_start = pos
            elif codon_type == 'stop' and last_start is not None:
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