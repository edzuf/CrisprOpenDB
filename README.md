# CrisprOpenDB

### Prerequisites

Before using this program, download the spacer database in the working directory (http://crispr.genome.ulaval.ca/dash/SpacersDB.fasta).

Then go to the `SpacersDB/` directory and download the sqlite file (http://crispr.genome.ulaval.ca/dash/CrisprOpenDB.sqlite). 

Once these steps are complete, you must go back to the previous directory to run the program.

### Running

To run the program, you must launch the `CrisperOpenDB_HostID.py` file in the working directory.
Here is an example of how to run the program using a *Salmonella* phage genome and a number of mismatches of 2:
```python
python CrisprOpenDB_HostID.py -i Salmonella_161.fasta -m 2
```
### Options

Alignment can be done using `blastn` or `fasta36`. If using Blast, please use `makeblastdb` before running. Here is the command line you should use when running `makeblastdb`:
```python
makeblastdb -in SpacersDB.fasta -dbtype nucl -out SpacersDB
```
