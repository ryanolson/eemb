from Queue import Queue, PriorityQueue
from threading import Thread
from decimal import Decimal
from copy import copy, deepcopy

import os
import re
import time
import itertools


class Job(object):
    def __init__(self, priority, description, nodes=1, ppn=1):
        self.priority = priority
        self.description = description
        self.nodes = nodes
        self.ppn = ppn
        return
    def type(self):
        return "Abstract Job"
    def cores(self):
        return self.nodes*self.ppn
    def __cmp__(self, other):
        return cmp(self.priority, other.priority)

class GamessCalculation(Job):
    def __init__(self, priority, description, nodes=1, ppn=1, qm=None):
        super(GamessCalculation, self).__init__(priority, description, nodes, ppn)
        if qm:
           self.qm = deepcopy(qm)
        return
    def type(self):
        return "GAMESS"
    def template_string(self):
        return """
 $CONTRL coord=unique SCFTYP=rhf runtyp=energy cctyp=ccsd(t) $END
 $contrl maxit=60 ispher=1 iskprp=1 icut=11 nprint=-5 $end
 $system parall=.t. mwords=20 memddi=400 timlim=999999 $end
 $BASIS GBASIS=acct $end
 $scf fdiff=.f. diis=.t. soscf=.f. dirscf=.t. $end
 $ccinp ikcut=1 ncore=0 maxcc=50 $end
 $DATA
 $CCSD(T) - {desc}
 C1
""".format(desc=self.description)
    def coordinates(self):
        s = ""
        for m in self.qm:
            if s != "":
               s += "\n"
            s += m.coordinates()
        return s
    def input_string(self):
        return self.template_string() + self.coordinates() + "\n $end\n\n"
    def update_results(self, results):
        """
 Results to extract from a GAMESS CCSD(T) calculation

      REFERENCE ENERGY:     -116.6000062236
        MBPT(2) ENERGY:     -117.0794036021   CORR.E=  -0.4793973784
        CCSD    ENERGY:     -117.1516954144   CORR.E=  -0.5516891908
        CCSD[T] ENERGY:     -117.1516954144   CORR.E=  -0.5516891908
        CCSD(T) ENERGY:     -117.1516954144   CORR.E=  -0.5516891908

 Timing Breakdown

 TOTAL WALL CLOCK TIME=        2.5 SECONDS, CPU UTILIZATION IS 100.00%
 TOTAL WALL CLOCK TIME=        2.6 SECONDS, CPU UTILIZATION IS 100.00%
 TOTAL WALL CLOCK TIME=        3.3 SECONDS, CPU UTILIZATION IS 100.00%
 TOTAL WALL CLOCK TIME=        3.4 SECONDS, CPU UTILIZATION IS 100.00%
 TOTAL WALL CLOCK TIME=      204.9 SECONDS, CPU UTILIZATION IS 100.00%
 TOTAL WALL CLOCK TIME=      257.2 SECONDS, CPU UTILIZATION IS 100.00%
 ITER WALL TIME: TOTAL=   704.0  TERMS IN AO BASIS=   680.5, MO BASIS=    23.4
 TOTAL WALL CLOCK TIME=     4011.9 SECONDS, CPU UTILIZATION IS 100.00%
 TOTAL WALL CLOCK TIME=     4114.6 SECONDS, CPU UTILIZATION IS 100.00%

"""
        pass  
    def run(self, threadID, nidlist=None, results=None):
        wb = re.compile(r'\s+')
        path = "data/{0}".format(wb.sub('-', self.description))
        mkdir(path)
        inp = "{0}/input.inp".format(path)
        out = "{0}/output.log".format(path)
        #
        # if the directory contains an output, update the results
        # and if the run does not need to be re-launch, return
        #
        if os.path.exists(out):
           return
           self.update_results(results)
           if results.get('status',None):
              return
        # create input file
        f = open(inp, "w")
        f.write(self.input_string())
        f.close()
        # launch job
        nidstr = ",".join(nidlist)
        cmd = "cd {path}; ln -f -s ../../rungms-xt; ./rungms-xt input.inp asis {cores} {ppn} \"{nid}\" {threadID}".format(path=path, cores=self.cores(), ppn=self.ppn, nid=nidstr, threadID=threadID)
        import subprocess
        subprocess.call( cmd, shell=True )
        # os.system( cmd )
        # import subprocess as sub
        # p = sub.Popen( cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True )
        # output, errors = p.communicate()
        # f = open(out, "w")
        # f.write( output )
        # f.close()


class AtomicCoordinate(object):
    def __init__(self, znuc, x, y, z, label=None):
        self.znuc = int(znuc)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.l = str(label)
        return
    def label(self):
        if self.l:
           return self.l
        return "RMO"
    def coordinate_string(self):
        return "{l}    {zn}.0   {x:15.10f} {y:15.10f} {z:15.10f}".format(l=self.label(), zn=self.znuc, x=self.x, y=self.y, z=self.z)


class Monomer(object):
    def __init__(self, description, atoms):
        self.description = description
        self.atoms = deepcopy(atoms)
        return
    def coordinates(self):
        s = ""
        for a in self.atoms:
            if s != "":
               s += "\n"
            s += a.coordinate_string()
        return s

def distanceBetween(m1, m2):
    x1 = y1 = z1 = 0.0
    x2 = y2 = z2 = 0.0
    for a in m1.atoms:
        x1 += a.znuc*a.x
        y1 += a.znuc*a.y
        z1 += a.znuc*a.z
    for a in m2.atoms:
        x2 += a.znuc*a.x
        y2 += a.znuc*a.y
        z2 += a.znuc*a.z
    x1 /= len(m1.atoms)
    y1 /= len(m1.atoms)
    z1 /= len(m1.atoms)
    x2 /= len(m2.atoms)
    y2 /= len(m2.atoms)
    z2 /= len(m2.atoms)
    r2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
    return r2


def mkdir(d):
    if not os.path.exists(d):
       os.makedirs(d)

def JobLauncher(threadID, q, nodeq):
    while True:
        job = q.get()
        nids = [ ]
        for i in range(job.nodes):
            nid = nodeq.get()
            nids.append(nid)
        print """
        {thread}: Launching new {type} job on {nodes} nodes using {ppn} cores per node [{cores} total cores]
                  Description: {desc}
                  Nidlist: {nids}
        """.format(thread=threadID,type=job.type(),nodes=job.nodes, ppn=job.ppn, cores=job.cores(), desc=job.description,nids=nids)
        job.run(threadID, nidlist=nids,results=results)
        q.task_done()
        for nid in nids:
            nodeq.put( nid )


def readCoordinateFile(coord_file):
    m = [ ] 
    wb = re.compile(r'\s+')
    with open(coord_file, "r") as f:
       for line in f:
           line = wb.sub(' ', line.strip()) 
           (label, znuc, x, y, z) = line.split(' ')
           m.append( AtomicCoordinate(int(Decimal(znuc)), x, y, z, label) )
           if len(m) == 3:
              monomers.append( Monomer("water {0}".format(len(monomers)+1), m) )
              m = [ ]

def generateNidList(nodes):
    cmd = "aprun -n {N} -N1 ./mpi_nidlist".format(N=nodes)
    os.system( cmd )

def parseNidList():
    f = open("nidlist", "r")
    for line in f:
        nid = line.strip()
        nidlist.append(nid)
        nodeq.put( nid )
    f.close()
    thread_count = len(nidlist)
    print thread_count

    

nidlist      = [ ]
monomers     = [ ]
results      = { }
q = PriorityQueue()
nodeq = Queue()

def initializeThreads(thread_count):
    for i in range(thread_count):
        worker = Thread(target=JobLauncher, args=(i,q,nodeq))
        worker.setDaemon(True)
        worker.start()

# move to unittest
# m1 = Monomer("monomer 1", [
#              AtomicCoordinate( 8.0, -1.2983607766,  -.5634207960, -4.0744440181, label="O" ),
#              AtomicCoordinate( 1.0, -1.3342246790, -1.3745953028, -3.5177484740, label="H" ),
#              AtomicCoordinate( 1.0, -1.3856343708,  -.9050814482, -4.9861409241, label="H" )
#      ])
# print m1.coordinates()

# move to unittest
# q.put( GamessCalculation( 4, "monomer") )
# q.put( GamessCalculation( 2, "trimer 1") )
# q.put( GamessCalculation( 2, "trimer 2") )
# q.put( GamessCalculation( 2, "trimer 3") )
# q.put( GamessCalculation( 2, "trimer 4") )
# q.put( GamessCalculation( 2, "trimer 5", nodes=5, ppn=8) )
# q.put( GamessCalculation( 3, "dimer") )
# q.join()

if __name__ == '__main__':
    import sys
    start = time.time()

    # input args
    argv = sys.argv
    node_count = int(Decimal(argv[1]))
    thread_count = int(Decimal(argv[2]))
    coord_file  = argv[3]

    # read input
    readCoordinateFile(coord_file)

    # generate and read nidlist
    generateNidList(node_count)
    parseNidList()

    # queue up trimers
    for trimer in itertools.combinations(monomers, 3):
        description = "{m1}-{m2}-{m3}".format(m1=trimer[0].description, m2=trimer[1].description, m3=trimer[2].description)
        td = 0
        td += distanceBetween(trimer[0], trimer[1])
        td += distanceBetween(trimer[0], trimer[2])
        td += distanceBetween(trimer[1], trimer[2])
        priority = int(td)
        calc = GamessCalculation( priority, description, nodes=4, ppn=8, qm=trimer )
        q.put( calc )
    print "{0} total calculations".format(q.qsize())

    # initialize threads & launch calculations
    initializeThreads(thread_count)
    q.join()

    finish = time.time()
    print "walltime = {time}".format(time=(finish-start))
    print "nodeq.qsize() = {0}".format(nodeq.qsize())

