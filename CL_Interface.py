from CrisprOpenDB.CrisprOpenDB_HostID import PhageHostFinder
import argparse

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file in FASTA format.", type=str, required=True)
    parser.add_argument("-m", "--mismatch", help="Number of mismatches. Value must be between 0 and 5. If not specified, default value is 1.", type=int, default=1)
    parser.add_argument("-a", "--aligner", help="Alignement tool to use. blast or fasta36", type=str, default="blast")
    parser.add_argument("-r", "--report", help="Show report of host identification.", action="store_true")
    parser.add_argument("-t", "--table", help="Show full result table in separate csv file.", action="store_true")
    parser.add_argument("-b", "--blastdb", help="Blast database to use for alignment.", type=str, default=None)
    parser.add_argument("-f", "--fastadb", help="Fasta database to use for alignment.", type=str, default=None)
    parser.add_argument("-u", "--unknown", help="Keep spacers with unknown genus for prediction. False if not specified.", action="store_true")
    args = parser.parse_args()

    if args.mismatch < 0 or args.mismatch > 5:
        parser.print_help()
        exit()
    
    if args.aligner not in ["blast", "fasta36"]:
        parser.print_help()
        exit()

    if args.unknown:
        print("**Warning**\nKeeping spacers with unknown genus (option -u, --unknown) can lead to incorrect prediction. If you choose to use this option, we stongly recommend that you use the report of host identification to keep track of the identification process and avoid any bias (option -r, --report).")

    if (args.blastdb != None and args.fastadb != None):
        print("Please use only one of the following options:\n-b, --blastdb\n-f, --fastadb")   
        exit()
    elif args.blastdb:
        phf = PhageHostFinder(args.blastdb, None)
    elif args.fastadb:
        phf = PhageHostFinder(None, args.fastadb)
    else:
        phf = PhageHostFinder()

    results = phf.identify(args.input, args.mismatch, args.aligner, args.report, args.table, args.unknown)
