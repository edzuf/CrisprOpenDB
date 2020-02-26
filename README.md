# CrisprOpenDB

### Prerequisites

Before using this program, go to the `CrisprOpenDB/` directory and download the spacer database (http://crispr.genome.ulaval.ca/dash/SpacersDB.fasta).

Then go the `CrisprOpenDB/CrisprOpenDB/` directory and dowload the sqlite file (http://crispr.genome.ulaval.ca/dash/CrisprOpenDB.sqlite). 

### Execution

To run the program, you must launch the `CrisperOpenDB_HostID.py` file in the `CrisprOpenDB` directory.
Here is an example of how run the program using a *Salmonella* phage genome and 2 mismatches:
```python
python CrisprOpenDB_HostID.py -i Salmonella_161.fasta -m 2
```
### Options

Alignment can be done using `blastn` or `fasta36`.
