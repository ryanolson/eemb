from monomer import readCoordinateFile

def eemb_efpfit(args):
    print "EFP Fitting Module"
    print args
    monomers = readCoordinateFile(args.xyzfile)

    if args.plot_fit_type == 'dft':
       analyzeDFTFit(args,monomers)
       pass

    if args.plot_fit_type == 'rhf':
       analyzeRHFFit(args,monomers)
       pass

def analyzeDFTFit(args,monomers):
    for m in monomers:
        if(args.show_qm_coords): print m.coordinates()
        if(args.show_efp_coords): print m.efp_dft_coordinates()
        val = diff(m.atoms, m.efp_dft)
        print m.description, val

def analyzeRHFFit(args,monomers):
    for m in monomers:
        if(args.show_qm_coords): print m.coordinates()
        if(args.show_efp_coords): print m.efp_dft_coordinates()
        val = diff(m.atoms, m.efp_rhf)
        print m.description, val

def diff(qm, efp):
    if len(qm) != len(efp):
       raise EEMBSizeMismatch
    val = 0
    for i in range(len(qm)):
       val += qm[i].distanceBetween(efp[i])
    return val
