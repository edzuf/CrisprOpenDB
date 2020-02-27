import os, sys, shlex, subprocess, io, base64, time
from Bio import SeqIO, Seq
from collections import Counter
from uuid import uuid4
import pandas as pd
import numpy as np
import argparse
from SpacersDB import CrisprOpenDB

class PhageHostFinder:
    def __init__(self):
        self._blast_database = "SpacersDB"
        self._fasta_database = "SpacersDB.fasta"
        self._alignement_results = None
        self._spacer_table= None


    def _run_fasta_36(self, fasta_file):
        command = "fasta36 -m 8 {} {}".format(fasta_file, self._database)
        print("Running command... {}".format(command))
        command = shlex.split(command)
        p = subprocess.Popen(command, stdout=subprocess.PIPE)
        p.wait()
        p = p.stdout.read()
        columns = ["Query", "SPACER_ID", "identity", "alignement_length", 
        "mismatch", "gap", "q_start", "q_end", "s_start", "s_end", "e_value", "score"]

        fasta36_out = io.StringIO(p.decode("utf-8"))
        fasta36_out.seek(0)
        fasta_result_table = pd.read_table(fasta36_out, names=columns)

        if len(fasta_result_table) == 0:
            print("No hits found. Sorry!")
            sys.exit()
        self._alignement_results = fasta_result_table
    
    def _run_blastn(self, fasta_file):
        command = "blastn -task 'blastn' -query {} -db {} -outfmt 6 ".format(fasta_file, self._blast_database)
        print("Running command... {}".format(command))
        command = shlex.split(command)
        p = subprocess.Popen(command, stdout=subprocess.PIPE)
        p.wait()
        p = p.stdout.read()

        columns = ["Query", "SPACER_ID", "identity", "alignement_length", 
        "mismatch", "gap", "q_start", "q_end", "s_start", "s_end", "e_value", "score"]

        blastn_out = io.StringIO(p.decode("utf-8"))
        blastn_out.seek(0)
        blastn_result_table = pd.read_table(blastn_out, names=columns)

        if len(blastn_result_table) == 0:
            print("No hits found. Sorry!")
            sys.exit()
        self._alignement_results = blastn_result_table
    
    def _load_spacer_table(self):
        t1 = time.time()
        db_explorer = CrisprOpenDB.CrisprOpenDB(os.path.join("SpacersDB", "CrisprOpenDB.sqlite"))
        df = pd.read_sql_query("select ST.SPACER_ID, ST.GENEBANK_ID, ORG.ORGANISM_NAME, ORG.SPECIES, ORG.GENUS, ORG.FAMILY, ORG.TORDER, ST.SPACER, ST.SPACER_LENGTH, SAL.COUNT_SPACER, ST.POSITION_INSIDE_LOCUS  \
            from ORGANISM ORG, SPACER_TABLE ST, SPACER_ARRAY_LENGTH SAL \
            where ST.GENEBANK_ID=ORG.GENEBANK_ID and ST.GENEBANK_ID=SAL.GENEBANK_ID and ST.NUMERO_LOCUS=SAL.NUMERO_LOCUS", db_explorer._connection)
        print("Database query took: {}".format(time.time()-t1))
        self._spacer_table = df


    @property
    def alignement_results(self):
        for query in set(np.array(self._alignement_results["Query"])):
           yield self._alignement_results[self._alignement_results["Query"] == query]
          
    def _findHost(self, alignement_table, n_mismatch):
        if self._spacer_table is None: #Lazy loading. Only load if needed (if there are alignement results).
            self._load_spacer_table()

        
        fasta_result_table = pd.merge(alignement_table, self._spacer_table, on="SPACER_ID", how="left")

        fasta_result_table["true_num_mismatch"] = (fasta_result_table["SPACER_LENGTH"] - fasta_result_table["alignement_length"]) + fasta_result_table["mismatch"]
        fasta_result_table  = fasta_result_table[fasta_result_table.true_num_mismatch <= n_mismatch]

        fasta_result_table["mean_position"] = np.array((fasta_result_table["q_start"] + fasta_result_table["q_end"]) / 2, dtype=int)

        #Criteria 1: If only one genus, it is the host.
        if len(set(fasta_result_table["GENUS"])) == 1:
            #print("Host is {}. Found using criteria #1".format(fasta_result_table["GENUS"][0]))
            return({"Query": fasta_result_table["Query"][0],
                "Host": fasta_result_table["GENUS"][0], 
                "Level": 1})
        
        #Criteria 2: Number of different position on phage genome (use mean(start, end) to check positions) (ex: MF153391)
        genus = np.array(fasta_result_table["GENUS"])
        position = np.array(fasta_result_table["mean_position"])
        sets_to_count = list(set([(j, position[i]) for i,j in enumerate(genus)]))
        counted_genus = Counter([i[0] for i in sets_to_count])

        most_commons_genus = counted_genus.most_common()
        print(most_commons_genus)
        if most_commons_genus[0][1] != most_commons_genus[1][1]: # If count is not equal, we found host
            #print("Host is {}. Found using criteria #2".format(most_commons_genus[0][0]))
            return({"Query": fasta_result_table["Query"][0],
                "Host": most_commons_genus[0][0], 
                "Level": 2})
        genuses_to_keep = []
        for i, g in enumerate(most_commons_genus):
            if most_commons_genus[i][1] == most_commons_genus[0][1]:
                genuses_to_keep.append(g[0])
        print("Genus to keep:" + str(genuses_to_keep))
        
        fasta_result_table = fasta_result_table[fasta_result_table["GENUS"].isin(genuses_to_keep)]

        #Criteria 3: If 2 is equal, find relative position most in 5' (MF158036)
        t1 = time.time()
        five_prime_relative_position = (np.array(fasta_result_table["POSITION_INSIDE_LOCUS"]) -1) / (np.array(fasta_result_table["COUNT_SPACER"]) - 1)
        fasta_result_table["five_prime_relative_position"] = five_prime_relative_position
        fasta_result_table.sort_values(by="five_prime_relative_position", inplace=True)
        fasta_result_table.reset_index(inplace=True, drop=True)
        #fasta_result_table.to_csv("TEST_CRITERIA3.csv")
        sub_section = fasta_result_table[fasta_result_table["five_prime_relative_position"] ==fasta_result_table["five_prime_relative_position"][0]]
        if len(sub_section==1):
            #print("Host is {}. Found using criteria #3".format(sub_section["GENUS"][0]))
            return({"Query": fasta_result_table["Query"][0],
                "Host": sub_section["GENUS"][0], 
                "Level": 3})
        else:
            fasta_result_table = sub_section

        #Criteria 4: Last common ancester (does not return a genus)
        print("Criteria 3: {}".format(time.time()-t1))

        if len(set(fasta_result_table["FAMILY"])) == 1:
            return({"Query": fasta_result_table["Query"][0],
                "Host": fasta_result_table["FAMILY"][0], 
                "Level": 4})
        elif len(set(fasta_result_table["TORDER"])) == 1:
            return({"Query": fasta_result_table["Query"][0],
                "Host": fasta_result_table["TORDER"][0], 
                "Level": 4})
        else:
            return({"Query": fasta_result_table["Query"][0],
                "Host": "UNKNOWN", 
                "Level": 4})


    def identify(self, fasta_file, n_mismatch, tool="blast"):
        if tool == "blast":
            self._run_blastn(fasta_file)
        elif tool == "fasta36":
            self._run_fasta_36(fasta_file)
        
        for table in self.alignement_results:
            hostIdentification = self._findHost(table, n_mismatch)
            print(hostIdentification)
        return(hostIdentification)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file in FASTA format.", required=True)
    parser.add_argument("-m", "--mismatch", help="Number of mismatches. Value must be between 0 and 5. If not specified, default value is 1.", type=int, default=1)
    parser.add_argument("-a", "--aligner", help="Alignement tool to use. blast or fasta36", type=str, default="blast")
    args = parser.parse_args()

    if args.mismatch < 0 or args.mismatch > 5:
        parser.print_help()
        exit()
    
    if args.aligner not in ["blast", "fasta36"]:
        parser.print_help()
        exit()

    phf = PhageHostFinder()
    results = phf.identify(args.input, args.mismatch, args.aligner)

    print(results)


