# bio-info-lab1
Current functionality:
  * Reading fasta format files and finding all longest start/stop codon pairs


A file to use is passed as a launch argument 
```bash
py find_orfs.py viruses/data/bacterial1.fasta
```
Multiple files can be processed by giving a folder path (non recursive)
```bash
py find_orfs.py viruses/data/
```

Currently there are two methods of finding codon pairs implemented
The alternative can be run by using --method argument
```bash
py find_orfs.py viruses/data/ --method alternative
```
The differances between the methods can be compared by using --compare  (comparison in a single file, folder cannot be passed)
```bash
py find_orfs.py viruses/data/bacterial1.fasta --compare
```
