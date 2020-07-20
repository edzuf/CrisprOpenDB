import os, sys, shlex, subprocess, io, base64, time
from Bio import SeqIO, Seq
from collections import Counter
from uuid import uuid4
import pandas as pd
import numpy as np
import argparse
from CrisprOpenDB.SpacersDB import CrisprOpenDB

class PhageHostFinder:
    def __init__(self, blast_db=os.path.join("CrisprOpenDB", "SpacersDB", "SpacersDB"), fasta_db=os.path.join("CrisprOpenDB", "SpacersDB", "SpacersDB.fasta")):
        self._blast_database = blast_db
        self._fasta_database = fasta_db
        self._alignement_results = None
        self._connection = None


    def _run_fasta_36(self, fasta_file, num_threads):
        if self._fasta_database is None: #TODO: add check if file exists
            print("Please specify a fasta database to perform analysis using fasta36.")
            sys.exit()
        elif not os.path.exists(self._fasta_database):
            print("Fasta database not found.")
            sys.exit()
        command = "fasta36 -T {} -m 8 {} {}".format(num_threads, fasta_file, self._fasta_database)
        #print("Running command... {}".format(command))
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
    
    def _run_blastn(self, fasta_file, num_threads):
        if self._blast_database is None: #TODO: add check if file exists
            print("Please specify a formated BLAST database to perform analysis using BLAST")
            sys.exit()

        list_file = self._blast_database.rsplit("/", 1)
        if len(list_file) == 1:
            directory = "."
        else:
            directory = list_file[0]
        f_exists = False
        for f in os.listdir(directory):
            if f.startswith(list_file[-1]):
                f_exists = True
                break
        if not f_exists:
            print("BLAST database not found.")
            sys.exit()
        
        command = "blastn -task 'blastn' -query {} -db {} -num_threads {} -outfmt 6".format(fasta_file, self._blast_database, num_threads)
        # print("Running command... {}".format(command))
        command = shlex.split(command)
        p = subprocess.Popen(command, stdout=subprocess.PIPE)
        p, err = p.communicate()
        columns = ["Query", "SPACER_ID", "identity", "alignement_length", 
        "mismatch", "gap", "q_start", "q_end", "s_start", "s_end", "e_value", "score"]

        blastn_out = io.StringIO(p.decode("utf-8"))
        blastn_out.seek(0)
        
        blastn_result_table = pd.read_csv(blastn_out, names=columns, sep="\t")

        blastn_out.close()

        if len(blastn_result_table) == 0:
            print("No hits found. Sorry!")
            sys.exit()
        self._alignement_results = blastn_result_table

    @property
    def alignement_results(self):
        for query in set(np.array(self._alignement_results["Query"])):
           yield self._alignement_results[self._alignement_results["Query"] == query]
          
    def _findHost(self, alignement_table, n_mismatch, report, table_to_file, keep_unknown):
        fasta_result_table = alignement_table.copy()

        #1. Select spacers id from the fasta result table
        spacers_id_from_blast = np.array(fasta_result_table.index)
        chunks = [spacers_id_from_blast[x:x+900] for x in range(0, len(spacers_id_from_blast), 900)]

        bd_df = pd.DataFrame()

        #2. Extract data from DB
        for element in chunks:
            temp_df = pd.read_sql_query(sql="select ST.SPACER_ID, ST.GENEBANK_ID, ORG.ORGANISM_NAME, ORG.SPECIES, ORG.GENUS, ORG.FAMILY, \
                ORG.TORDER, ST.SPACER, ST.SPACER_LENGTH, SAL.COUNT_SPACER, ST.POSITION_INSIDE_LOCUS \
                from ORGANISM ORG, SPACER_TABLE ST, SPACER_ARRAY_LENGTH SAL \
                where ST.GENEBANK_ID=ORG.GENEBANK_ID \
                and ST.GENEBANK_ID=SAL.GENEBANK_ID \
                and ST.NUMERO_LOCUS=SAL.NUMERO_LOCUS \
                and ST.SPACER_ID in ({seq})".format(seq=','.join(['?']*len(element))), 
                params=element, 
                con=self._connection._connection)
            if bd_df.empty:
                bd_df=temp_df
            else:
                bd_df=pd.concat([bd_df, temp_df])
        
        #3. Merge.
        fasta_result_table = fasta_result_table.merge(bd_df, on="SPACER_ID", how="left")

        print("\n====================")
        print("Query: {}".format(fasta_result_table["Query"].iloc[0]))
        print("====================")

        if not keep_unknown:
            fasta_result_table = fasta_result_table[fasta_result_table["GENUS"] != "Unknown"]
    
        fasta_result_table.loc[:,"true_num_mismatch"] = (fasta_result_table.loc[:,"SPACER_LENGTH"] - fasta_result_table.loc[:,"alignement_length"]) + fasta_result_table.loc[:,"mismatch"]

        if table_to_file:
            fasta_result_table.to_csv("{}.csv".format(fasta_result_table["Query"].iloc[0]))

        fasta_result_table = fasta_result_table[fasta_result_table.gap <= 0]
        fasta_result_table  = fasta_result_table[fasta_result_table.true_num_mismatch <= n_mismatch]

        if fasta_result_table.empty:
            if report:
                print("Empty result table. Found no match between query and spacer database that tolerates {} mismatch(es).".format(n_mismatch))
            return("No hits found. Sorry!")

        fasta_result_table["mean_position"] = np.array((fasta_result_table["q_start"] + fasta_result_table["q_end"]) / 2, dtype=int)

        #Criteria 1: If only one genus, it is the host.
        if report:
            print("\n**Criteria 1: If only one genus, it is the host**\n")

        if len(set(fasta_result_table["GENUS"])) == 1:
            if report:
                print("Host is {}. Found using criteria #1".format(fasta_result_table["GENUS"].iloc[0]))
                print("Spacer(s) on which prediction is based:")
                for i in range(len(fasta_result_table["SPACER"])):
                    print("{}; organism: {}; start position on phage genome: {}; end position on phage genome: {}".format(fasta_result_table["SPACER"].iloc[i], 
                    fasta_result_table["ORGANISM_NAME"].iloc[i], fasta_result_table["q_start"].iloc[i], fasta_result_table["q_end"].iloc[i]))
                
            return({"Query": fasta_result_table["Query"].iloc[0],
                "Host": fasta_result_table["GENUS"].iloc[0],  
                "Level": 1})
        if report:
            print("Multiple possible genuses:")
            for i in set(fasta_result_table["GENUS"]):
                print(i)

        #Criteria 2: Number of different position on phage genome (use mean(start, end) to check positions) (ex: MF153391)
        if report:
            print("\n**Criteria 2: Number of different positions on phage genome**\n")

        genus = np.array(fasta_result_table["GENUS"])
        position = np.array(fasta_result_table["mean_position"])
        sets_to_count = list(set([(j, position[i]) for i,j in enumerate(genus)]))
        counted_genus = Counter([i[0] for i in sets_to_count])

        most_commons_genus = counted_genus.most_common()
        
        if report:
            for i in range(len(most_commons_genus)):
                print("Genus: {}; different positions: {}".format(most_commons_genus[i][0], most_commons_genus[i][1]))
        if most_commons_genus[0][1] != most_commons_genus[1][1]: # If count is not equal, we found host
            if report:
                print("Host is {}. Found using criteria #2".format(most_commons_genus[0][0]))
                print("Spacer(s) on which prediction is based:")
                for i in range(len(fasta_result_table["SPACER"])):
                    if fasta_result_table["GENUS"].iloc[i] == most_commons_genus[0][0]:
                        print("{}; organism: {}; start position on phage genome: {}; end position on phage genome: {}".format(fasta_result_table["SPACER"].iloc[i],
                        fasta_result_table["ORGANISM_NAME"].iloc[i], fasta_result_table["q_start"].iloc[i], fasta_result_table["q_end"].iloc[i]))
            return({"Query": fasta_result_table["Query"].iloc[0],
                "Host": most_commons_genus[0][0], 
                "Level": 2})
        genuses_to_keep = []
        for i, g in enumerate(most_commons_genus):
            if most_commons_genus[i][1] == most_commons_genus[0][1]:
                genuses_to_keep.append(g[0])
        
        if report:
            print("Genuses to keep:" + str(genuses_to_keep))
        fasta_result_table = fasta_result_table[fasta_result_table["GENUS"].isin(genuses_to_keep)]

        #Criteria 3: If 2 is equal, find relative position most in 5' (MF158036)
        if report:
            print("\n**Criteria 3: If number of positions on phage genome is equal for multiple genuses, find relative position most in 5'**\n")
    
        five_prime_relative_position = (np.array(fasta_result_table["POSITION_INSIDE_LOCUS"], dtype=float) -1) / (np.array(fasta_result_table["COUNT_SPACER"]) - 1)
        fasta_result_table["five_prime_relative_position"] = five_prime_relative_position
        fasta_result_table.sort_values(by="five_prime_relative_position", inplace=True)
        fasta_result_table.reset_index(inplace=True, drop=True)
        
        if report:
            print("5' relative positions for all remaining potential hosts:")
            print(fasta_result_table[["GENUS", "five_prime_relative_position"]])
        sub_section = fasta_result_table[fasta_result_table["five_prime_relative_position"] == fasta_result_table["five_prime_relative_position"][0]]
        
        if len(sub_section["GENUS"].unique().tolist()) == 1:
            if report:
                print("Host is {}. Found using criteria #3".format(sub_section["GENUS"].iloc[0]))
                print("Spacer(s) on which prediction is based:")
                for i in range(len(sub_section["SPACER"])):
                    print("{}; organism: {}; start position on phage genome: {}; end position on phage genome: {}".format(sub_section["SPACER"].iloc[i], 
                    sub_section["ORGANISM_NAME"].iloc[i], sub_section["q_start"].iloc[i], sub_section["q_end"].iloc[i]))
            return({"Query": fasta_result_table["Query"].iloc[0],
                "Host": sub_section["GENUS"].iloc[0], 
                "Level": 3})
        else:
            fasta_result_table = sub_section
            if report:
                print("Potential hosts with the same relative position (most in 5'): ")
                for genus in set(fasta_result_table["GENUS"]):
                    print(genus)

        #Criteria 4: Last common ancester (does not return a genus)
        
        if report:
            print("\n**Criteria 4: Last common ancester (does not return a genus)**\n")
            print(fasta_result_table[["GENUS", "FAMILY", "TORDER"]])

        if len(set(fasta_result_table["FAMILY"])) == 1:
            if report:
                print("Family is common to all remaining potential hosts.")
                print("Host is {}. Found using criteria #4".format(fasta_result_table["FAMILY"].iloc[0]))
                print("Spacer(s) on which prediction is based:")
                for i in range(len(fasta_result_table["SPACER"])):
                    print("{}; organism: {}, start position on phage genome: {}; end position on phage genome: {}".format(fasta_result_table["SPACER"].iloc[i], 
                    fasta_result_table["ORGANISM_NAME"].iloc[i], fasta_result_table["q_start"].iloc[i], fasta_result_table["q_end"].iloc[i]))
            return({"Query": fasta_result_table["Query"].iloc[0],
                "Host": fasta_result_table["FAMILY"].iloc[0], 
                "Level": 4})
        elif len(set(fasta_result_table["TORDER"])) == 1:
            if report:
                print("Order is common to all remaining potential hosts.")
                print("Host is {}. Found using criteria #4".format(fasta_result_table["TORDER"].iloc[0]))
                print("Spacer(s) on which prediction is based:")
                for i in range(len(fasta_result_table["SPACER"])):
                    print("{}; organism: {}, start position on phage genome: {}; end position on phage genome: {}".format(fasta_result_table["SPACER"].iloc[i], 
                    fasta_result_table["ORGANISM_NAME"].iloc[i], fasta_result_table["q_start"].iloc[i], fasta_result_table["q_end"].iloc[i]))
            return({"Query": fasta_result_table["Query"].iloc[0],
                "Host": fasta_result_table["TORDER"].iloc[0], 
                "Level": 4})
        else:
            if report:
                print("Unable to find last common ancester.")
                print("Spacer(s) on which prediction is based:")
                for i in range(len(fasta_result_table["SPACER"])):
                    print("{}; organism: {}, start position on phage genome: {}; end position on phage genome: {}".format(fasta_result_table["SPACER"].iloc[i], 
                    fasta_result_table["ORGANISM_NAME"].iloc[i], fasta_result_table["q_start"].iloc[i], fasta_result_table["q_end"].iloc[i]))
            return({"Query": fasta_result_table["Query"].iloc[0],
                "Host": "UNKNOWN", 
                "Level": 4})


    def identify(self, fasta_file, n_mismatch, tool="blast", report=False, table_to_file=False, keep_unknown=False, num_threads=1):
        if tool == "blast":
            self._run_blastn(fasta_file, num_threads)
        elif tool == "fasta36":
            self._run_fasta_36(fasta_file, num_threads)
        
        if len(self._alignement_results) == 0:
            return("No hits found. Sorry!")
        
        df =  self._alignement_results.set_index("SPACER_ID")
        self._alignement_results = df

        if self._connection is None:
            self._connection = CrisprOpenDB.CrisprOpenDB(os.path.join("CrisprOpenDB", "SpacersDB", "CrisprOpenDB.sqlite"))
        
        for table in self.alignement_results:
            hostIdentification = self._findHost(table, n_mismatch, report, table_to_file, keep_unknown)
            print("\n{}".format(hostIdentification))

        return(hostIdentification)



