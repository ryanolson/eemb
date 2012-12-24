from monomer import readCoordinateFile

def eemb_analyze(args):
    print "EEMB Analyzer"
    monomers = readCoordinateFile(args.xyzfile)

    if args.process_data:
       pass

    if args.with_summary:
       pass


