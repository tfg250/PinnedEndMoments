from subprocess import CREATE_NEW_CONSOLE
import subprocess
from subprocess import STDOUT
import copy
import os
import pandas as pd
import math
import re
import shelve
import plotly.offline as py
import plotly.graph_objs as go
import pprint
import csv
pp = pprint.PrettyPrinter(indent=4)



def recovershelve(OutputFile, skips=[], only=False):
    '''
    function to recover a shelve contents
    outputfile is the name of the shelf (without date & time)
    skips are the contents to skip 
    only can be specified to collect only certain variables
    '''
    # find the most recent shelve file...
    path = os.getcwd()
    files = os.listdir(path)
    shelves = []

    # add date regex
    pattern1 = re.compile(
        "[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]")

    for f in files:
        if f[-3:] == 'shv' and OutputFile in f and bool(re.search(pattern1, f)):
            shelves.append(f)
    if shelves == []:
        if os.path.isfile(OutputFile+".shv"):
            shelves.append(OutputFile+".shv")

    shelffile = sorted(shelves)[-1]
    #print "recovering "+shelffile
    RETURNS = {}

    my_shelf = shelve.open(shelffile, 'r')
    for key in sorted(my_shelf):
        if not only:
            if key not in ['bridgemiserstart', 'timepoint1', 'timepoint0', 'timeelapsed', 'Logger', 'Tee', 'color_cycle', 'FILE', 'overwriteshelveratinginputs', 'SKIPSHELVES', 'DEBUG_MODE', 'DEBUGFILE', 'Loggertimestring'] and key not in skips:
                try:
                    RETURNS[key] = my_shelf[key]
                except:
                    printred("Cannot unshelve key: "+key)
        else:
            if key in only:
                try:
                    RETURNS[key] = my_shelf[key]
                    #print key
                except:
                    printred("Cannot unshelve key: "+key)
    my_shelf.close()
    return RETURNS

def GetLogFile(OutputFile):
    path = os.getcwd()
    files = os.listdir(path)
    logfiles = []
    # add date regex
    pattern1 = re.compile(
        "[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]")

    for f in files:
        if f[-3:] == 'log' and OutputFile in f and bool(re.search(pattern1, f)) and 'logger' not in f:
            logfiles.append(f)
    if logfiles == []:
        if os.path.isfile(OutputFile+".shv"):
            logfiles.append(OutputFile+".shv")
    if len(logfiles) == 0:
        return 'None'
    else:
        logfile = sorted(logfiles)[-1]
        return logfile


from txdot_inputs import *



MODELDATA = {}
recoverkeys = ['GIRDERLINEELEMENTS', 'NODES', 'BEAMS', 'BCnodes', 'G_listofelements',
               'NODES_OUT', 'G_listofys', 'COMP', 'NONCOMP', 'REACTIONS', 'DISPLACEMENTS', 'NODEMAP']



#read csv of runs alredy processed
csvin = 'txdot_results.csv'
input_file = csv.DictReader(open(csvin))
RESULTS= {}
for row in input_file:
    model = row["Model"]
    RESULTS[model]={k:row[k] for k in row.keys() if k!="Model"}


success = 0
finished = []
for inp in [inp for inp in TX if inp not in RESULTS]:
    #print inp,
    # check that it completed successfully
    logfile = GetLogFile(inp)
    if os.path.isfile(logfile):
        with open(logfile, 'r') as f:
            if "ELAPSED SOLUTION TIME" in f.readlines()[-2]:
                MODELDATA[inp] = recovershelve(inp, [], recoverkeys)
                if 'REACTIONS' not in MODELDATA[inp]:
                    print "ERROR IN ",inp," MISSING REACTIONS"
                if 'COMP' not in MODELDATA[inp]:
                    print "ERROR IN ",inp," MISSING COMP"
                if 'DISPLACEMENTS' not in MODELDATA[inp]:
                    print "ERROR IN ",inp," MISSING DISPLACEMENTS"
                if 'REACTIONS' in MODELDATA[inp] and 'COMP' in MODELDATA[inp] and 'DISPLACEMENTS' in MODELDATA[inp]:
                    finished.append(inp)
                    success += 1
                    print '\t',inp,
                    print "Success"

            if not success:
                print '\t',inp,
                print "Failed"
    else:
        print '\t',inp,
        print "not yet"


print success




#build dataframe of results
if True:
    for inp in finished:
        RESULTS[inp]= {'width':int(inp[3:5]),'span':TX[inp][1],'girderspacing':TX[inp][2],'numberofgirders':TX[inp][3],'girders':TX[inp][4],'numberofdiaphragms':TX[inp][5],'Y':TX[inp][6]}
        maxXreaction = 0
        minXreaction = 0
        for n in MODELDATA[inp]['REACTIONS']:
            cases = [c for c in MODELDATA[inp]['REACTIONS'][n].keys() if c not in ['DL2', 'DL1', 'DL0', 'SIDL']]
            maxXreaction = max(maxXreaction,max([MODELDATA[inp]['REACTIONS'][n][c][0] for c in cases]))
            minXreaction = min(minXreaction,min([MODELDATA[inp]['REACTIONS'][n][c][0] for c in cases]))
        RESULTS[inp]['maxXreaction'] = maxXreaction
        RESULTS[inp]['minXreaction'] = minXreaction
        maxmoment = 0
        mindeflection = 0
        for yg in sorted(MODELDATA[inp]['GIRDERLINEELEMENTS'].keys()):
            for tup in MODELDATA[inp]['GIRDERLINEELEMENTS'][yg]:
                g = 'g'+str(tup[2])
                section = "SC_"+"%04d" % tup[0]
                cases = [c for c in MODELDATA[inp]['COMP'][g][section].keys() if c not in ['DL2', 'DL1', 'DL0', 'SIDL']]
                maxmoment = max(maxmoment, float(max([MODELDATA[inp]['COMP'][g][section][c][4] for c in cases])))
                ni = MODELDATA[inp]['NODEMAP'][MODELDATA[inp]['BEAMS'][tup[0]][2]]
                nj = MODELDATA[inp]['NODEMAP'][MODELDATA[inp]['BEAMS'][tup[0]][3]]
                mindeflection = min(mindeflection,float(min(min([MODELDATA[inp]['DISPLACEMENTS'][ni][c][2] for c in cases]),min([MODELDATA[inp]['DISPLACEMENTS'][nj][c][2] for c in cases]))))
        RESULTS[inp]['maxmoment'] = maxmoment
        RESULTS[inp]['mindeflection'] = mindeflection
        maxendmoment = 0
        minendmoment = 0
        for yg in sorted(MODELDATA[inp]['GIRDERLINEELEMENTS'].keys()):
            for tup in [MODELDATA[inp]['GIRDERLINEELEMENTS'][yg][0],MODELDATA[inp]['GIRDERLINEELEMENTS'][yg][-1]]:
                g = 'g'+str(tup[2])
                section = "SC_"+"%04d" % tup[0]
                cases = [c for c in MODELDATA[inp]['COMP'][g][section].keys() if c not in ['DL2', 'DL1', 'DL0', 'SIDL']]
                maxendmoment = max(maxendmoment, (float(max([MODELDATA[inp]['COMP'][g][section][c][4] for c in cases]))))
                minendmoment = min(minendmoment, (float(min([MODELDATA[inp]['COMP'][g][section][c][4] for c in cases]))))
        RESULTS[inp]['maxendmoment'] = maxendmoment
        RESULTS[inp]['minendmoment'] = minendmoment

resDF = pd.DataFrame.from_dict(RESULTS,orient='index')
resDF.to_csv('txdot_results.csv')





# Add calc of girder A & I
# other params needed in csv?


