# bio-info-lab1
Current functionality - reading fasta format files and finding all longest start/stop codon pairs, filtering them by length, converting into proteins, analysing their frequency and finaly producing two distance matrices in Phylip format, one for codones andone for dicodones.

A file to use is passed as a launch argument, this outputs all proteins found in the genome
```bash
py find_orfs.py viruses/data/bacterial1.fasta
```
Multiple files can be processed by giving a folder path (non recursive), use this method to output the distance matrices
```bash
py find_orfs.py viruses/data/
```

Currently there are two methods of finding codon pairs implemented
The alternative can be run by using --method argument
```bash
py find_orfs.py viruses/data/ --method alternative
```
Both methods find the same results
