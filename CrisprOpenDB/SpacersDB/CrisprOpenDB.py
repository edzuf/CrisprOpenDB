import sqlite3
import sys, os, re, glob, time
from Bio import Entrez

Entrez.email="pier-luc.plante.1@ulaval.ca"

class CrisprOpenDB:
    def __init__(self, db_file = "CrisprOpenDB.sqlite"):
        self.db_file = db_file
        self._connection = sqlite3.connect(db_file)
        self._cursor = self._connection.cursor()

    def create_tables(self, overwrite=False):
        if overwrite:
            self._cursor.execute("drop table if exists SPACER_TABLE")
            self._cursor.execute("drop table if exists ORGANISM")
            self._connection.commit()
        
        self._cursor.execute("create table if not exists ORGANISM (\
            GENEBANK_ID text, \
            ORGANISM_NAME text, \
            PRIMARY KEY (GENEBANK_ID), \
            SPECIES text, \
            GENUS text, \
            FAMILY text, \
            TORDER text)")
        self._cursor.execute("create table if not exists SPACER_TABLE (\
            SPACER_ID text, \
            SPACER text, \
            START integer, \
            END integer, \
            SPACER_LENGTH integer, \
            STRAND text(1), \
            POSITION_INSIDE_LOCUS integer, \
            NUMERO_LOCUS integer, \
            GENEBANK_ID text, \
            PRIMARY KEY (SPACER_ID), \
            FOREIGN KEY (GENEBANK_ID) REFERENCES ORGANISM(GENEBANK_ID))")
        self._connection.commit()
        return

    def add_new_organism(self, genebank, organism):
        self._cursor.execute("insert into ORGANISM (GENEBANK_ID, ORGANISM_NAME) values (?, ?)",
                            (genebank, organism))
        self._connection.commit()
    
    def insert_information(self, spacer, start, end, spacer_length, strand, genebank, pos_locus, num_locus):
        spacer_id = genebank + "_" + str(num_locus) + "_" + str(pos_locus) #We could use UUID instead.
        self._cursor.execute("insert into SPACER_TABLE values (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            (spacer_id, spacer, start, end, spacer_length, strand, pos_locus, num_locus, genebank))
        #self._connection.commit()
        return spacer_id

    def NCBIEntrez(self, accession):
        time.sleep(0.4)
        Entrez.email = "moira.dion.1@ulaval.ca"
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        next(handle)
        species = " ".join((handle.readline().split())[1:])
        return species
    
    def ExtractSpacers(self, inFile):
        with open(inFile, 'r') as fl:
            locus_count = 0
            accession = ""
            for line in fl:
                array = line.strip().split('\t')
                if array[2] == 'repeat_region':
                    locus_count += 1
                    spacer_count = 0
                    if array[0] != accession:
                        try:
                            accession = array[0]
                            #print("Asking NCBI for %s" %(accession))
                            species = self.NCBIEntrez(accession)
                        except Exception as e:
                            print(e)
                            yield([-1, accession])
                        try:
                            self.add_new_organism(accession, species)
                        except:
                            pass
                    else:
                        accession = accession
                        species = species
                if array[2] == 'binding_site':
                    #unique_spacer_id = 
                    spacer_count += 1
                    start_position= int(array[3])
                    end_position= int(array[4])
                    size= int(array[5])
                    orientation = array[6]
                    sArray = array[8].split(';')
                    for it in sArray:
                        if re.match('Note=', it):
                            sequence = it[5:]
                    yield([locus_count, spacer_count, start_position, end_position, size, orientation, sequence, accession])
            self._connection.commit()
    
    def fill_tables(self, path):
        print(path)
        file_list = glob.glob(os.path.join(path, "*.gff"))
        file_count=0
        t1 = time.time()
        for f in file_list:
            print(f)
            for spacer in self.ExtractSpacers(f):
                if spacer[0] == -1:
                    with open("DB_Creation_errors", 'a') as fo:
                        fo.write(f + "\t" + spacer[1]+"\n")
                    break
                self.insert_information(spacer[6], spacer[2], spacer[3], spacer[4], spacer[5], spacer[7], spacer[1], spacer[0])
            file_count += 1
            if not file_count%10:
                print("Done %f4%% in %f2" %(float(file_count)/len(file_list) * 100, time.time() - t1 ))
        return

    def update_table(self, file_in):
        """
        Add a collection of sequence to the database.
        file_in is a text file containing a list of gff file to process.
        """
        for line in open(file_in, 'r'):
            f = line.split()[0]
            print(f)
            for spacer in self.ExtractSpacers(f):
                if spacer[0] == -1:
                    with open("DB_update_errors", 'a') as fo:
                        fo.write(f + "\t" + spacer[1]+"\n")
                    break
                spacer_id = spacer[7] + "_" + str(spacer[0]) + "_" + str(spacer[1])
                self._cursor.execute("select 1 from SPACER_TABLE where SPACER_ID=? and \
                SPACER=? and \
                START=? and \
                END=? and \
                SPACER_LENGTH=? and \
                STRAND=? and \
                GENEBANK_ID=?", 
                (spacer_id, spacer[6], spacer[2], spacer[3], spacer[4], spacer[5], spacer[7]))
                if len(self._cursor.fetchall()) == 0:
                    try:
                        self.insert_information(spacer[6], spacer[2], spacer[3], spacer[4], spacer[5], spacer[7], spacer[1], spacer[0])
                    except Exception as e:
                        print(e)
        return

    def get_spacers_from_sequence(self, genebank):
        self._cursor.execute("select SPACER from SPACER_TABLE where GENEBANK_ID = ?", (genebank,)) #La "," du deuxieme parametre est necessaire
        query_resuls = [i[0] for i in self._cursor.fetchall()]
        return(query_resuls)

    def get_spacers_from_sequences_iterator(self, genebanks):
        for genebank in genebanks:
            self._cursor.execute("select SPACER from SPACER_TABLE where GENEBANK_ID = ?", (genebank,))
            for sequence in self._cursor:
                yield(genebank, sequence)

    def create_complete_fasta_file(self, filename):
        with open(filename, 'w') as fi:
            self._cursor.execute("select SPACER_ID, SPACER from SPACER_TABLE")
            for spacer_entry in self._cursor:
                fi.write(">%s\n" %spacer_entry[0])
                fi.write(spacer_entry[1]+"\n")
    
    def create_subset_fasta_file(self, organism_keyword):
        self._cursor.execute("select SPACER_ID, SPACER from SPACER_TABLE \
            where GENEBANK_ID in \
                (select GENEBANK_ID from ORGANISM where ORGANISM_NAME like ?)", ("%" + organism_keyword + "%",))
        query_resuls = [i[0] for i in self._cursor.fetchall()]
        return query_resuls

        
    def get_genebanks_from_organism(self, organism_keyword):
        self._cursor.execute("select GENEBANK_ID from ORGANISM where ORGANISM_NAME like ?)", ("%" + organism_keyword + "%",))
        query_resuls = [i[0] for i in self._cursor.fetchall()]
        return query_resuls 

    def close(self):
        self._connection.commit()
        self._connection.close()

    def count_number_of_spacers(self, distinct=False):
        if distinct:
            self._cursor.execute("select count(distinct SPACER) from SPACER_TABLE")
        else:
            self._cursor.execute("select count(SPACER_ID) from SPACER_TABLE ")
        return(self._cursor.fetchall()[0])

    def count_number_of_organisme(self):
        self._cursor.execute("select count(GENEBANK_ID) from ORGANISM")
        return(self._cursor.fetchall()[0])

    def create_view_for_app(self):
        print("Creating view...")
        self._cursor.execute("create view ORGANISM_SPACER as \
             select SPACER_TABLE.GENEBANK_ID, ORGANISM_NAME, SPACER, SPACER_LENGTH  \
                from ORGANISM, SPACER_TABLE \
                where SPACER_TABLE.GENEBANK_ID=ORGANISM.GENEBANK_ID")
        self._connection.commit()

    def fill_taxonomy_columns(self):
        self._cursor.execute("select distinct ORGANISM_NAME from ORGANISM")
        query_results = [i[0] for i in self._cursor.fetchall()]
        print("{} items to query".format(len(query_results)))
        counter = 0
        for org_id in query_results:
            time.sleep(1)
            try:
                handle = Entrez.esearch(db="Taxonomy", term=org_id)
                record = Entrez.read(handle)
                if len(record["IdList"]) == 1:
                    time.sleep(1)
                    print("{} record number is {}".format(org_id, record["IdList"][0]))
                    handle = Entrez.efetch(db="Taxonomy", id=record["IdList"][0], retmode="xml")
                    records = Entrez.read(handle)
                    species = ""
                    genus = ""
                    order = ""
                    suborder = ""
                    family = ""
                    for items in records[0]["LineageEx"]:
                        if items["Rank"] == "species":
                            species = items["ScientificName"]
                        elif items["Rank"] == "genus":
                            genus = items["ScientificName"]
                        elif items["Rank"] == "order":
                            order = items["ScientificName"]
                        elif items["Rank"] == "suborder":
                            suborder = items["ScientificName"]
                        elif items["Rank"] == "family":
                            family = items["ScientificName"]
                    print("Order: {}; Suborder: {}; family: {}, genus: {}, species: {}".format(order, suborder, family, genus, species))
                    self._cursor.execute("update ORGANISM set SPECIES = ?, \
                        GENUS = ?, \
                        FAMILY = ?, \
                        SUBORDER = ?, \
                        TORDER = ? \
                        where ORGANISM_NAME = ?", (species, genus, family, suborder, order, org_id))
                    counter += 1
                    if counter%10 == 0:
                        print("Commiting to DB")
                        self._connection.commit()
            except:
                with open("Error_taxonomy.txt", "a") as ferror:
                    print("Error with {}".format(org_id))
                    ferror.write("{}\n".format(org_id))
        self._connection.commit()
            

                
if __name__ == "__main__":
    test = CrisprOpenDB()
    test.fill_taxonomy_columns()


