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
pp = pprint.PrettyPrinter(indent=4)


def ComputeSectionProps(ttf,btf,tbf,bbf,tw,hw):
    """
    Function to compute composite section properties
    Args:
        ttf (float): thickness of top flange
        btf (float): width of top flange
        tbf (float): thickness of bottom flange
        bbf (float): width of bottom flange
        d (float): depth of steel section
        tw (float): thickness of web
        haunch (float): haunch
        bef (float): effective width (transformed into equivalent width of steel)
        ts (float): thickness of slab
    Returns:
        ybar (float): distance from bottom of steel to centroid (L)
        I (float): moment of inertia of the section (L**4)
        St (float): section modulus to the top of steel (L*3)
        Sb (float): section modulus to the bottom of steel (L**3)
    """
    d=hw+ttf-tbf
    A=tbf*bbf+ttf*btf+tw*hw
    ybar=(tbf*bbf*(tbf/2.) + tw*hw*(tbf+hw/2.) + ttf*btf*(tbf+hw+ttf/2.))/A
    return ybar

def ScrubLogFile(logfilename, masssum=False):
    """
    Function to determine the success of a run
    Args:
        logfilename (str): name of the log file
    Returns:
        numberofnotes (str): number of notes
        numberofwarnings (str): number of warnings
        numberoferrors (str): number of errors
        solution (str): completed?
        TimesOfLessThanMaxReaction (list): times at which the full truck isnt on the span
    """
    LogRead = open(logfilename, 'r')
    LineRead = LogRead.readline()
    solution = 'fail'
    numberofnotes = 'X'
    numberofwarnings = 'X'
    numberoferrors = 'X'
    FZList = []
    TimeList = []
    try:
        while LineRead:
            if LineRead.split(":")[0].strip() == "*Number of Notes":
                numberofnotes = LineRead.split(":")[1].strip()
            if LineRead.split(":")[0].strip() == "*Number of Warnings":
                numberofwarnings = LineRead.split(":")[1].strip()
            if LineRead.split(":")[0].strip() == "*Number of Errors":
                numberoferrors = LineRead.split(":")[1].strip()
            if "*Solution terminated" in LineRead:
                solution = 'Terminated'
            if "*Solution completed" in LineRead:
                solution = 'Completed'
            if "DIRECT SUMMATION OF NODE REACTION FORCES" in LineRead:
                #print LineRead.split()
                try:
                    timestep = int(LineRead.split()[8])
                except IndexError:
                    timestep = 0
                LineRead = LogRead.readline()
                LineRead = LogRead.readline()
                FZ = float(LineRead.split()[3])
                FZList.append(FZ)
                TimeList.append(timestep)
            LineRead = LogRead.readline()
        LogRead.close()
    except:
        pass
    #print logfilename,TimesOfLessThanMaxReaction
    if masssum:
        return numberofnotes, numberofwarnings, numberoferrors, solution, maxreaction
    else:
        return numberofnotes, numberofwarnings, numberoferrors, solution


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


# -------------------------------------------------------------------------------
# Write input file
# -------------------------------------------------------------------------------
def WriteInput(INPUTS):
    # now build the input file(s)
    print "Writing "+INPUTS['ModelName']+".inp"
    inpfilename = INPUTS['ModelName']+".inp"
    inpfile = open(inpfilename, 'w')

    inpfile.write("*OUTPUTFILE,"+INPUTS['ModelName']+",\n")
    inpfile.write("*UNITS,"+",".join([str(j) for j in INPUTS['Units']])+",\n")
    inpfile.write("*OPTIONS,"+",".join([str(j)
                                        for j in INPUTS['Options']])+",\n")
    inpfile.write("*DECKMESH,"+str(INPUTS['deckmesh'])+",\n")
    for m in sorted(INPUTS['Mats'].keys()):
        inpfile.write("*MATERIAL,"+str(int(m))+"," +
                      ",".join([str(j) for j in INPUTS['Mats'][m]])+",\n")

    for bp in sorted(INPUTS['Props'].keys()):
        inpfile.write("*BEAMPROP,"+str(int(bp))+"," +
                      ",".join([str(j) for j in INPUTS['Props'][bp]])+",\n")

    # DECK-------------------------------------------------------------------
    d1W = (INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] -
                                      1) + 2 * INPUTS['Overhang']) - 2 * max(INPUTS['SidewalkWidth'],INPUTS['ParapetWidth'])
    d1L = INPUTS['SpanLength']

    # assumed D1 is at z=0
    DL0 = 0
    DL1 = 0
    DL2 = 0
    DW = 0
    KF = 0
    inpfile.write("*DECK,1,"+str(d1L)+','+str(d1W)+',0,'+str(INPUTS['DeckMat'])+','+str(
        INPUTS['DeckThick'])+',0,'+str(-d1W/2.)+','+','.join([str(j) for j in [DL0, DL1, DL2, DW, KF]])+"\n")

    d2W = INPUTS['SidewalkWidth']
    d2L = INPUTS['SpanLength']


    y0=-(INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] -1) + 2 * INPUTS['Overhang'])/2.
    y1=y0+min(INPUTS['SidewalkWidth'],INPUTS['ParapetWidth'])
    y2=y0+max(INPUTS['SidewalkWidth'],INPUTS['ParapetWidth'])

    y5=(INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] -1) + 2 * INPUTS['Overhang'])/2.
    y3=y5-max(INPUTS['SidewalkWidth'],INPUTS['ParapetWidth'])
    y4=y5-min(INPUTS['SidewalkWidth'],INPUTS['ParapetWidth'])
    #print y0,y1,y2,y3,y4,y5
    #print (y1-y0)>0, (y2-y1)>0, (y4-y3)>0, (y5-y4)>0
    if INPUTS['SidewalkThick'] > 0 and INPUTS['SidewalkStiffness']:
        t_sidewalk = INPUTS['DeckThick']+INPUTS['SidewalkThick']
        Dw=0
        z = (INPUTS['DeckThick'] + INPUTS['SidewalkThick']) / 2. - INPUTS['DeckThick'] / 2.
    elif INPUTS['SidewalkThick'] > 0 and not INPUTS['SidewalkStiffness']:
        t_sidewalk = INPUTS['DeckThick']
        DWsw = INPUTS['SidewalkThick'] * INPUTS['Mats'][INPUTS['DeckMat']][2]
        z=0
    else:
        t_sidewalk = INPUTS['DeckThick']
        DWsw = 0
        z = 0
    if INPUTS['ParapetThick'] > 0 and not INPUTS['ParapetStiffness']:
        DWp = INPUTS['ParapetThick'] * INPUTS['Mats'][INPUTS['DeckMat']][2]
    else:
        DWp=0

    # ADD DECK 2,3,4,5 as needed
    if (y1-y0)>0:
        #2
        inpfile.write("*DECK,2,"+str(d2L)+','+str(y1-y0)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(
            t_sidewalk)+',0,'+str(y0)+','+','.join([str(j) for j in [DL0, DL1, DL2, DWsw+DWp, KF]])+"\n")
    if (y2-y1)>0:
        #3
        inpfile.write("*DECK,3,"+str(d2L)+','+str(y2-y1)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(
            t_sidewalk)+',0,'+str(y1)+','+','.join([str(j) for j in [DL0, DL1, DL2, DWsw, KF]])+"\n")
    if (y4-y3)>0:
        #4
        inpfile.write("*DECK,4,"+str(d2L)+','+str(y4-y3)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(
            t_sidewalk)+',0,'+str(y3)+','+','.join([str(j) for j in [DL0, DL1, DL2, DWsw, KF]])+"\n")
    if (y5-y4)>0:
        #5
        inpfile.write("*DECK,5,"+str(d2L)+','+str(y5-y4)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(
            t_sidewalk)+',0,'+str(y4)+','+','.join([str(j) for j in [DL0, DL1, DL2, DWsw+DWp, KF]])+"\n")

    # GIRDERS-------------------------------------------------------------------
    gx = 0
    y0 = -(INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] - 1))/2.
    BCs = []
    BCe = []
    for g in range(INPUTS['NumberGirders']):
        # *GIRDER,No,X1,X2,Y,Haunch,Prop,Mat,DL0,DL1,DL2,DW,Comp,Deck,KF,Tfubl,Bfubl,TSS,Cope
        if g in [0, INPUTS['NumberGirders']-1]:
            # exterior
            for x in range(int((len(INPUTS['Girder_E'])-1)/2.)):
                gx += 1
                xstart = INPUTS['Girder_E'][2*x]*d1L
                prop = INPUTS['Girder_E'][2*x+1][0]
                haunch = INPUTS['Girder_E'][2*x+1][1]
                xend = INPUTS['Girder_E'][2*x+2]*d1L
                inpfile.write("*GIRDER,"+str(gx)+','+str(xstart)+"," +
                              str(xend)+','+str(y0)+','+str(haunch)+','+str(prop))
                # Girder_common = [haunch, MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
                inpfile.write(
                    ','+','.join([str(j) for j in INPUTS['Girder_Ecommon']])+"\n")
                if x == 0:
                    BCs.append(gx)
                if x == int((len(INPUTS['Girder_E'])-1)/2.)-1:
                    BCe.append(gx)
        else:
            # interior
            for x in range(int((len(INPUTS['Girder_I'])-1)/2.)):
                gx += 1
                xstart = INPUTS['Girder_I'][2*x]*d1L
                prop = INPUTS['Girder_I'][2*x+1][0]
                haunch = INPUTS['Girder_I'][2*x+1][1]
                xend = INPUTS['Girder_I'][2*x+2]*d1L
                inpfile.write("*GIRDER,"+str(gx)+','+str(xstart)+"," +
                              str(xend)+','+str(y0)+','+str(haunch)+','+str(prop))
                # Girder_common = [haunch, MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
                inpfile.write(
                    ','+','.join([str(j) for j in INPUTS['Girder_Icommon']])+"\n")
                if x == 0:
                    BCs.append(gx)
                if x == int((len(INPUTS['Girder_E'])-1)/2.)-1:
                    BCe.append(gx)
        y0 += INPUTS['GirderSpacing']

    # compute deck inputs
    if INPUTS['ParapetWidth'] > 0 and INPUTS['ParapetThick'] > 0 and INPUTS['ParapetStiffness']:
        # make a new beam prop
        inpfile.write("*BEAMPROP,"+str(int(999))+',0,SolidTrapezoid,'+str(
            INPUTS['ParapetWidth'])+','+str(INPUTS['ParapetWidth'])+','+str(INPUTS['ParapetThick'])+'\n')
        # now the girders
        yp = ((INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] -
                                          1) + 2 * INPUTS['Overhang']) - INPUTS['ParapetWidth'])/2.
        haunch = -INPUTS['DeckThick']-INPUTS['ParapetThick']
        if INPUTS['SidewalkThick'] > 0 and INPUTS['SidewalkWidth'] > 0 and INPUTS['SidewalkStiffness']:
            haunch -= INPUTS['SidewalkThick']
        inpfile.write("*GIRDER,998,0,"+str(d1L)+','+str(-yp)+',' +
                      str(haunch)+',999,'+str(INPUTS['DeckMat'])+',0,0,0,0,B,1,3\n')
        inpfile.write("*GIRDER,999,0,"+str(d1L)+','+str(yp)+',' +
                      str(haunch)+',999,'+str(INPUTS['DeckMat']) + ',0,0,0,0,B,1,3\n')
        # *GIRDER, Girder Number, Xstart, Xend, Ycg, haunch, BeamPropNo, MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP (C;NC;NA), corresponding deck no, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE

    # *BC_GIRDER, Girder Number, X, Z, DOF1, DOF2, DOF3, DOF4, DOF5, DOF6 ,,,,,,,,,
    for gs in BCs:
        inpfile.write("*BC_GIRDER,"+str(gs)+',START,' +
                      str(INPUTS['BCZ'])+','+','.join([str(j) for j in INPUTS['BCS']])+"\n")
    for ge in BCe:
        inpfile.write("*BC_GIRDER,"+str(ge)+',END,' +
                      str(INPUTS['BCZ'])+','+','.join([str(j) for j in INPUTS['BCE']])+"\n")

    # DIAPHRAGMS-------------------------------------------------------------------
    # *DIAPHRAGM,No, zi, zm, zj,Mat,Prop,X, type,DL0,DL1,DL2,DW,KF,skip,,,,
    dno = 0
    for d in INPUTS['Diaphragm']:
        dno += 1
        inpfile.write('*DIAPHRAGM,'+str(dno)+','+str(d[3])+','+str(d[3])+','+str(
            d[3])+','+str(d[1])+','+str(d[2])+','+str(d[0]*d1L)+',')
        inpfile.write(','.join([str(j)
                                for j in INPUTS['Diaphragm_Common']])+'\n')

    if len(INPUTS['PointLoad']) > 0:
        lc=100
        
        for ptcase in INPUTS['PointLoad']:
            lc+=1
            inpfile.write('*STATICLOADCASE,'+str(lc)+',PointLoad'+str(INPUTS['PointLoad'].index(ptcase))+'\n')
            #print lc
            #pprint.pprint(ptcase)
            for ptload in ptcase:
                inpfile.write('*NODEFORCE,'+str(ptload[0]*d1L)+','+str(ptload[1]*d1W - d1W/2.)+',0,'+str(lc)+','+str(ptload[2])+','+str(ptload[3])+','+str(ptload[4])+',0,0,0\n')

    if len(INPUTS['Vehicles']) > 0:
        # *LANE, Lane Number, shape, curvature, x1, x2, y, z, width
        # Fit lanes between sidewalks
        if d1W <= 17*12:
            numlanes = 1
            if d1W <= 144:
                # one lane
                width = d1W
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,' +
                              str(d1L)+',0,0,'+str(width-20.)+'\n')
            else:
                width = 144
                # -
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
                # 0
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,' +
                              str(d1L)+',0,0,'+str(width)+'\n')
                # +
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANESETS, 1\n')
                inpfile.write('*LANESETS, 2\n')
                inpfile.write('*LANESETS, 3\n')

        elif d1W <= 308:
            numlanes = 2
            width = (d1W-48.)/2.

            # 0
            inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L) +
                          ','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
            inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L) +
                          ','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')
        else:
            numlanes = math.floor(d1W/144.)
            width = 144.-24.
            if numlanes == 2:
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-d1W/2.+3*width/2.+48.)+',0,'+str(width)+'\n')
                # 0
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-width/2.-24)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,4,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(width/2.+24)+',0,'+str(width)+'\n')
                # +
                inpfile.write('*LANE,5,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(d1W/2.-3*width/2.-48.)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,6,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANESETS, 1,2\n')
                inpfile.write('*LANESETS, 3,4\n')
                inpfile.write('*LANESETS, 5,6\n')

            elif numlanes == 3:
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(-d1W/2.+3*width/2.+48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-d1W/2.+5*width/2.+96)+',0,'+str(width)+'\n')
                # 0
                inpfile.write('*LANE,4,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-2*width/2.-24)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,5,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(0)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,6,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(2*width/2.+24)+',0,'+str(width)+'\n')
                # +
                inpfile.write('*LANE,7,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(d1W/2.-5*width/2.-96)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,8,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(d1W/2.-3*width/2.-48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,9,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')

                inpfile.write('*LANESETS, 1,2,3\n')
                inpfile.write('*LANESETS, 4,5,6\n')
                inpfile.write('*LANESETS, 7,8,9\n')
            else:
                if numlanes > 4:
                    print "limiting lanes to 4"
                    numlanes == 4
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(-d1W/2.+3*width/2.+48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(-d1W/2.+5*width/2.+96)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,4,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(-d1W/2.+7*width/2.+144)+',0,'+str(width)+'\n')
                # 0
                inpfile.write('*LANE,5,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-3*width/2.-72)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,6,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(-width/2.-24)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,7,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(width/2.+24)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,8,STRAIGHT,FLAT,0,'+str(d1L) +
                              ','+str(3*width/2.+72)+',0,'+str(width)+'\n')
                # +
                inpfile.write('*LANE,9,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(d1W/2.-7*width/2.-144)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,10,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(d1W/2.-5*width/2.-96)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,11,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(d1W/2.-3*width/2.-48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,12,STRAIGHT,FLAT,0,'+str(d1L)+',' +
                              str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')

                inpfile.write('*LANESETS, 1,2,3,4\n')
                inpfile.write('*LANESETS, 5,6,7,8\n')
                inpfile.write('*LANESETS, 9,10,11,12\n')

        for v in INPUTS['Vehicles']:
            inpfile.write('*VEHICLE,20,' +
                          str(v[0])+','+str(v[1])+','+str(v[2])+'\n')
        inpfile.write('*MPF,'+','.join([str(j) for j in INPUTS['MPF']])+'\n')
        inpfile.write('*IMPACTFACTOR,'+str(INPUTS['IM'])+'\n')
        inpfile.write('*LOADFACTOR,' +
                      ','.join([str(j) for j in INPUTS['LF']])+'\n')
        inpfile.write('*PHI,'+','.join([str(j) for j in INPUTS['PHI']])+'\n')

    inpfile.write('*PLOTLYMESH\n')
    #inpfile.write('*PLOTLYSOLID\n')
    #inpfile.write('*MULTI,6,6\n')

    inpfile.close()


# -------------------------------------------------------------------------------
# Run Miser
# -------------------------------------------------------------------------------
def RunMiser(INPUTS, runMODE, extractMODE, rateMODE, skipshelves=False):

    FileName = INPUTS['ModelName']+'.inp'
    try:
        if skipshelves:
            a = subprocess.check_call(["python32", 'bridgemiser_016.py', FileName, str(
                runMODE), str(extractMODE), str(rateMODE), 'on', 'on'], stderr=STDOUT)
        else:
            a = subprocess.check_call(["python32", 'bridgemiser_016.py', FileName, str(
                runMODE), str(extractMODE), str(rateMODE)], stderr=STDOUT)
    except Exception as e:
        #raise e
        return False
    return True


def MomentPlot(MODELDATA, INPUTS, run, filename, extract='static', norm=False):
    tableau10_256 = []
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120), (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150), (148, 103, 189), (197, 176, 213),
                 (140, 86, 75), (196, 156, 148), (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199), (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    for i in range(len(tableau20)/2):  # integer operation
        r, g, b = tableau20[i*2]
        #color_cycle2.append((r / 255., g / 255., b / 255.))
        tableau10_256.append((r, g, b))

    PlotTraces = []
    maxz = 0
    maxy = 0
    miny = 0
    maxx = 0
    minx = 0
    for inp in sorted(run):
        first = True
        cc = sorted(INPUTS.keys()).index(inp) % 10
        maxZinp = 0
        if norm:
            for yg in sorted(MODELDATA[inp]['GIRDERLINEELEMENTS'].keys()):
                for tup in MODELDATA[inp]['GIRDERLINEELEMENTS'][yg]:
                    g = 'g'+str(tup[2])
                    section = "SC_"+"%04d" % tup[0]
                    #pp.pprint(MODELDATA[inp]['COMP'].keys())
                    #pp.pprint(MODELDATA[inp]['COMP'][g][section].keys())
                    if extract=='static':
                        cases = [c for c in MODELDATA[inp]['COMP'][g][section].keys() if c not in ['DL2', 'DL1', 'DL0', 'SIDL']]
                        maxZinp = max(maxZinp, (float(max([MODELDATA[inp]['COMP'][g][section][c][4] for c in cases]))))
                    else:
                        maxZinp = max(maxZinp, (float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4])))


        for yg in sorted(MODELDATA[inp]['GIRDERLINEELEMENTS'].keys()):
            x = []
            y = []
            zmin = []
            zmax = []
            hovermax = []
            hovermin = []
            for tup in MODELDATA[inp]['GIRDERLINEELEMENTS'][yg]:
                # list of (bee,xi,g,xj)
                x.append(float(tup[1]))
                y.append(float(yg))
                g = 'g'+str(tup[2])
                section = "SC_"+"%04d" % tup[0]
                x.append(float(tup[3]))
                y.append(float(yg))
                if extract=='static':
                    cases = [c for c in MODELDATA[inp]['COMP'][g][section].keys() if c not in ['DL2', 'DL1', 'DL0', 'SIDL']]
                    zi1=float(max([MODELDATA[inp]['COMP'][g][section][c][4] for c in cases]))
                    zi2=float(min([MODELDATA[inp]['COMP'][g][section][c][4] for c in cases]))
                    zj1=float(max([MODELDATA[inp]['COMP'][g][section][c][4+6] for c in cases]))
                    zj2=float(min([MODELDATA[inp]['COMP'][g][section][c][4+6] for c in cases]))
                else:
                    cases = [6]
                    zi1=float(max([MODELDATA[inp]['COMP'][g][section][c]['MAX'][4] for c in cases]))
                    zi2=float(min([MODELDATA[inp]['COMP'][g][section][c]['MIN'][4] for c in cases]))
                    zj1=float(max([MODELDATA[inp]['COMP'][g][section][c]['MAX'][4+6] for c in cases]))
                    zj2=float(min([MODELDATA[inp]['COMP'][g][section][c]['MIN'][4+6] for c in cases]))
                if norm:
                    #pp.pprint(MODELDATA[inp]['COMP'][g][section].keys())
                    zmax.append(zi1/maxZinp)
                    zmin.append(zi2/maxZinp)
                    zmax.append(zj1/maxZinp)
                    zmin.append(zj2/maxZinp)
                    hovermax.append(zi1)
                    hovermin.append(zi2)
                    hovermax.append(zj1)
                    hovermin.append(zj2)
                else:
                    zmax.append(zi1)
                    zmin.append(zi2)
                    zmax.append(zj1)
                    zmin.append(zj2)
                #compoutputterms = ['Axial','HShear','VShear','Torsion','HMoment','VMoment']
            maxz = max(maxz, max(zmax))
            maxy = max(maxy, max(y))
            miny = min(miny, min(y))
            maxx = max(maxx, max(x))
            minx = min(minx, min(x))
            if norm:
                if first:
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmax,
                                           hoverinfo="all", text=hovermax, legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly'))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmin, hoverinfo="all",
                                           text=hovermin, legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly', showlegend=False))
                    first = False
                else:
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmax, hoverinfo="all",
                                           text=hovermax, legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly', showlegend=False))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmin, hoverinfo="all",
                                           text=hovermin, legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly', showlegend=False))
            else:
                if first:
                    PlotTraces.append(dict(type='scatter3d', line=dict(
                        color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmax, hoverinfo="all", legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly'))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmin,
                                           hoverinfo="all", legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly', showlegend=False))
                    first = False
                else:
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmax,
                                           hoverinfo="all", legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly', showlegend=False))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=zmin,
                                           hoverinfo="all", legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly', showlegend=False))

    noaxis = dict(showbackground=False, showline=False, zeroline=False,
                  showgrid=False, showticklabels=False, title='')
    layout = dict(
        title="Moment Summary", showlegend=True,
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(x=(maxx-minx)/(maxy-miny), y=1, z=1),
            camera=dict(eye=dict(x=0, y=-0.5, z=0)),
            xaxis=noaxis,
            yaxis=noaxis,
            zaxis=noaxis))
    fig = dict(data=PlotTraces, layout=layout)
    py.plot(fig, filename=filename+'.html', auto_open=False, show_link=False)


def DeflectionPlot(MODELDATA, INPUTS, run, filename, extract='static', SF=10000.):
    tableau10_256 = []
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120), (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150), (148, 103, 189), (197, 176, 213),
                 (140, 86, 75), (196, 156, 148), (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199), (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
    for i in range(len(tableau20)/2):  # integer operation
        r, g, b = tableau20[i*2]
        #color_cycle2.append((r / 255., g / 255., b / 255.))
        tableau10_256.append((r, g, b))

    if extract=='static':
        pattern1 = re.compile("PointLoad[0-9]")
    else:
        pattern1 = re.compile("_C06_L[0-9][0-9]")

    PlotTraces = []
    maxz = 0
    maxy = 0
    miny = 0
    maxx = 0
    minx = 0
    for inp in sorted(run):
        cases = []
        #print inp
        first = True
        cc = sorted(INPUTS.keys()).index(inp) % 10
        tmpnode = MODELDATA[inp]['DISPLACEMENTS'].keys()[0]
        for c in MODELDATA[inp]['DISPLACEMENTS'][tmpnode].keys():
            if bool(re.search(pattern1, str(c))):
                cases.append(c)
        case = sorted(cases)[0]
        #case = 6
        for yg in sorted(MODELDATA[inp]['GIRDERLINEELEMENTS'].keys()):
            x = []
            y = []
            z = []
            hover = []
            for tup in MODELDATA[inp]['GIRDERLINEELEMENTS'][yg]:
                # list of (bee,xi,g,xj)
                # x.append(float(tup[1]))
                # y.append(float(yg))
                g = 'g'+str(tup[2])
                section = "SC_"+"%04d" % tup[0]
                # x.append(float(tup[3]))
                # y.append(float(yg))
                #print MODELDATA[inp]['BEAMS'][tup[0]]
                #print MODELDATA[inp]['NODEMAP'][MODELDATA[inp]['BEAMS'][tup[0]][2]]
                #print MODELDATA[inp]['NODEMAP'][MODELDATA[inp]['BEAMS'][tup[0]][3]]
                ni = MODELDATA[inp]['NODEMAP'][MODELDATA[inp]
                                               ['BEAMS'][tup[0]][2]]
                nj = MODELDATA[inp]['NODEMAP'][MODELDATA[inp]
                                               ['BEAMS'][tup[0]][3]]
                #print ni in MODELDATA[inp]['DISPLACEMENTS']
                #print nj in MODELDATA[inp]['DISPLACEMENTS']
                #print ni in MODELDATA[inp]['NODES_OUT']
                #print nj in MODELDATA[inp]['NODES_OUT']
                #print MODELDATA[inp]['BEAMS'][tup[0]][2] in MODELDATA[inp]['NODES_OUT']
                #print MODELDATA[inp]['BEAMS'][tup[0]][3] in MODELDATA[inp]['NODES_OUT']
                #print MODELDATA[inp]['DISPLACEMENTS'][ni].keys()
                #print case
                #print case in MODELDATA[inp]['DISPLACEMENTS'][ni].keys()

                # + SF * float(MODELDATA[inp]['DISPLACEMENTS'][ni][6]['MIN'][0])
                x.append(float(MODELDATA[inp]['NODES_OUT']
                               [MODELDATA[inp]['BEAMS'][tup[0]][2]][0]))
                # + SF * float(MODELDATA[inp]['DISPLACEMENTS'][ni][6]['MIN'][1])
                y.append(float(MODELDATA[inp]['NODES_OUT']
                               [MODELDATA[inp]['BEAMS'][tup[0]][2]][1]))
                if extract=='static':
                    z.append(float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][2]][2]) + SF * float(MODELDATA[inp]['DISPLACEMENTS'][ni][case][2]))
                    hover.append('x: '+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][2]][0])+"<br>"+"y: "+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][2]][1])+"<br>" + "z: "+"%.2f" % float(MODELDATA[inp]['DISPLACEMENTS'][ni][case][2]))
                else:
                    z.append(float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][2]][2]) + SF * float(MODELDATA[inp]['DISPLACEMENTS'][ni][case]['MIN'][2]))
                    hover.append('x: '+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][2]][0])+"<br>"+"y: "+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][2]][1])+"<br>" + "z: "+"%.2f" % float(MODELDATA[inp]['DISPLACEMENTS'][ni][case]['MIN'][2]))
                # + SF * float(MODELDATA[inp]['DISPLACEMENTS'][ni][6]['MIN'][0])
                x.append(float(MODELDATA[inp]['NODES_OUT']
                               [MODELDATA[inp]['BEAMS'][tup[0]][3]][0]))
                # + SF * float(MODELDATA[inp]['DISPLACEMENTS'][ni][6]['MIN'][1])
                y.append(float(MODELDATA[inp]['NODES_OUT']
                               [MODELDATA[inp]['BEAMS'][tup[0]][3]][1]))
                if extract=='static':
                    z.append(float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][3]][2]) + SF * float(MODELDATA[inp]['DISPLACEMENTS'][nj][case][2]))
                    hover.append('x: '+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][3]][0])+"<br>"+"y: "+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][3]][1])+"<br>" + "z: "+"%.2f" % float(MODELDATA[inp]['DISPLACEMENTS'][nj][case][2]))
                else:
                    z.append(float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][3]][2]) + SF * float(MODELDATA[inp]['DISPLACEMENTS'][nj][case]['MIN'][2]))
                    hover.append('x: '+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][3]][0])+"<br>"+"y: "+"%.2f" % float(MODELDATA[inp]['NODES_OUT'][MODELDATA[inp]['BEAMS'][tup[0]][3]][1])+"<br>" + "z: "+"%.2f" % float(MODELDATA[inp]['DISPLACEMENTS'][nj][case]['MIN'][2]))

            maxz = max(maxz, -min(z))
            maxy = max(maxy, max(y))
            miny = min(miny, min(y))
            maxx = max(maxx, max(x))
            minx = min(minx, min(x))
            if first:
                PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=z,
                                       hoverinfo="name+text", text=hover, legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly'))
                first = False
            else:
                PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`, width=3), mode='lines', x=x, y=y, z=z, hoverinfo="name+text",
                                       text=hover, legendgroup=inp, name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'], visible='legendonly', showlegend=False))

    noaxis = dict(showbackground=False, showline=False, zeroline=False,
                  showgrid=False, showticklabels=False, title='')
    layout = dict(
        title="Deflection Summary", showlegend=True,
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(x=(maxx-minx)/(maxy-miny), y=1, z=1),
            camera=dict(eye=dict(x=0, y=-0.5, z=0)),
            xaxis=noaxis,
            yaxis=noaxis,
            zaxis=noaxis))
    fig = dict(data=PlotTraces, layout=layout)
    py.plot(fig, filename=filename+'.html', auto_open=False, show_link=False)


# -------------------------------------------------------------------------------
# Import Results
# -------------------------------------------------------------------------------
# def ReadResults(INPUTS):
#    compdemandsfile = INPUTS['ModelName']+'_comp_sectiondemands.csv'
#    reactionfile = INPUTS['ModelName']+'_nodalreactions.csv'
#    noncompdemandsfile = INPUTS['ModelName']+'_nodalreactions.csv'
#    displacementfile = INPUTS['ModelName']+"_nodaldisplacements.csv"
#    reactions = pd.read_csv(reactionfile)
#    displacements = pd.read_csv(displacementfile)
#    if os.path.isfile(compdemandsfile):
#        demands = pd.read_csv(compdemandsfile)
#    elif os.path.isfile(noncompdemandsfile):
#        demands = pd.read_csv(noncompdemandsfile)
#    return reactions, demands, displacements


# -------------------------------------------------------------------------------
# INPUTS
# -------------------------------------------------------------------------------
if False:
    #MBE_A1
    ModelName = 'MBEA1_P'
    gz = -24.0
    Units = ['in', 'lb', 'Btu', 'psi', 'lbf', 'F', 386.089]
    NumberGirders = 4
    GirderSpacing = 88.
    SidewalkWidth = 0
    ParapetWidth = 0
    SpanLength = 780.
    Overhang = 11
    DeckThick = 7.25
    DeckMat = 2
    SidewalkThick = 0.
    ParapetThick = 0.
    ParapetMat = 2
    SidewalkStiffness = True  # thickened shell
    ParapetStiffness = True  # makes a 'girder'
    # LoadLocation = [x,y,fx,fy,fz] #as % of length, width, 0.5,0.5 is center.
    #PointLoad = [[0.5,0.5,0,0,-1000.]]
    PointLoad = []
    Vehicles = [['HS20_14FT', 'DESIGN', 0.0]]
    MPF = [1.2, 1, 0.85, 0.65]
    IM = 1.33
    LF = [1.25,   0.9,     1.5,    0.65,         1.75,         1.35,      1.45,
          1.5,    1.0,    1.0,       1.3,      1.0,         1.3,      1.0]
    PHI = [1, 1]
    # x as % of span length, p as prop integer
    #          [x0,(prop,haunch),x1,(prop,haunch),x2,(prop,haunch),x3]
    Girder_E = [0, (1, 0), 0.2, (2, 0), 0.8, (1, 0), 1.0]
    Girder_I = [0, (1, 0), 0.2, (3, 0), 0.8, (1, 0), 1.0]
    # Girder_common = [MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
    Girder_Icommon = [1, 0, 0, 0, 0, 'C', 1, 1, 0, 195, 195, 0]
    Girder_Ecommon = [1, 0, 0, 0, 0, 'C', 1, 1, 0, 195, 195, 0]
    GirderType = 'C'
    #           [x,mat,prop,z]
    Diaphragm = [[0, 1, 4, gz], [0.25, 1, 4, gz], [
        0.5, 1, 4, gz], [0.75, 1, 4, gz], [1, 1, 4, gz]]
    Diaphragm_Common = ['straight', 0, 0, 0, 0, 1]
    Props = {}
    Props[1] = [0, 'IBeam', 11.5, 11.5, 33.1, 0.855, 0.855, 0.58]
    Props[2] = [0, 'IBeam', 11.0835, 11.5, 33.85, 1.48, 0.855, 0.58]
    Props[3] = [0, 'IBeam', 11.038, 11.5, 33.725, 1.605, 0.855, 0.58]
    Props[4] = [0, 'Lipchannel', 3.95, 18, 0, 0.625, 0.45, 0]

    # Props[1]=[1,0,Null,37.8712,6607.75,217.235,6.73246,40.75,20.375]
    # Props[2]=[2,0,Null,44.5148,8277.46,276.798,15.5185,40.75,20.375]
    # Props[3]=[3,0,Null,45.6822,8414.14,288.743,18.3722,40.75,20.375]
    # Props[4]=[4,0,Null,12.475,549.035,15.6833,1.10832,18,9]
    Mats = {}
    #Mats[no] = [E, nu, rho, strength]
    Mats[1] = [29.0e6, 0.29, 490./12.**3, 36000.]
    Mats[2] = [57000.*(4000.)**0.5, 0.29, 150./12.**3, 4000.]
    deckmesh = 10.
    skew = 0.
    # F for free, R for restrained, float for spring stiffness
    BCS = ['R', 'R', 'R', 'F', 'F', 'F']  # start of each girder
    BCE = ['F', 'F', 'R', 'F', 'F', 'F']  # end of each girder
    BCZ = 0.
    # guidance is truss elements, so release moments
    Options = ['diaphragmreleases', True,
               'reactions', True, 'displacements', True,'sectioncheckingfigure',False]



if False:
    #BCW57901
    #MBE_A1
    ModelName = 'BCW57901_P'
    gz = -30.745
    Units = ['in', 'lb', 'Btu', 'psi', 'lbf', 'F', 386.089]
    NumberGirders = 11
    GirderSpacing = 93.6
    SidewalkWidth = 0
    ParapetWidth = 12.
    SpanLength = 1284.0
    Overhang = 48
    DeckThick = 7.5
    DeckMat = 2
    SidewalkThick = 0.
    ParapetThick = 48.
    ParapetMat = 2
    SidewalkStiffness = True  # thickened shell
    ParapetStiffness = False  # makes a 'girder'
    # LoadLocation = [x,y,fx,fy,fz] #as % of length, width, 0.5,0.5 is center.
    #PointLoad = [[0.5,0.5,0,0,-1000.]]
    PointLoad = []
    Vehicles = [['HS20_14FT', 'DESIGN', 0.0]]
    MPF = [1.2, 1, 0.85, 0.65]
    IM = 1.33
    LF = [1.25,   0.9,     1.5,    0.65,         1.75,         1.35,      1.45,
          1.5,    1.0,    1.0,       1.3,      1.0,         1.3,      1.0]
    PHI = [1, 1]
    # x as % of span length, p as prop integer
    #          [x0,(prop,haunch),x1,(prop,haunch),x2,(prop,haunch),x3]
    Girder_E = [0, (1, 0), 234./1284., (2, 0), (1284.-234.)/1284., (1, 0), 1.0]
    Girder_I = [0, (1, 0), 234./1284., (2, 0), (1284.-234.)/1284., (1, 0), 1.0]
    # Girder_common = [MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
    Girder_Icommon = [1, 0, 0, 0, 0, 'C', 1, 1, 0, 1284/5., 1284/5., 0]
    Girder_Ecommon = [1, 0, 0, 0, 0, 'C', 1, 1, 0, 1284/5., 1284/5., 0]
    GirderType = 'C'
    #           [x,mat,prop,z]
    Diaphragm = [[0, 1, 5, gz-11.9], [0, 1, 6, gz+18.3], [0.2, 1, 4, gz], [0.4, 1, 4, gz], [0.6, 1, 4, gz], [0.8, 1, 4, gz], [1.0, 1, 5, gz-11.9], [1.0, 1, 6, gz+18.3]]
    Diaphragm_Common = ['straight', 0, 0, 0, 0, 1]
    Props = {}
    Props[1] = [0, 'IBeam',16,16,53.75,0.875,0.875,0.375]
    Props[2] = [0, 'IBeam', 16,16,54.5,1.625,0.875,0.375]
    Props[4] = [0, 'AISC', '13th','MC18x42.7']
    Props[5] = [0, 'AISC', '13th','C12x20.7']
    Props[6] = [0, 'AISC', '13th','C15x33.9']

    Mats = {}
    #Mats[no] = [E, nu, rho, strength]
    Mats[1] = [29.0e6, 0.29, 490./12.**3, 36000.]
    Mats[2] = [57000.*(3000.)**0.5, 0.29, 150./12.**3, 3000.]
    deckmesh = 10.
    skew = 0.
    # F for free, R for restrained, float for spring stiffness
    BCS = ['R', 'R', 'R', 'F', 'F', 'F']  # start of each girder
    BCE = ['F', 'F', 'R', 'F', 'F', 'F']  # end of each girder
    BCZ = 0.
    # guidance is truss elements, so release moments
    Options = ['diaphragmreleases', True,
               'reactions', True, 'displacements', True]
#                     0,  1,     2,3,                    4,5,    6]
#TX['TX_30_120_R01']=[2,120,     7,5,[40,277]             ,6,49.62]
#TX['TX_24_30_P01'] =[2, 30,7.3333,4,[1,12,1.25,12,0.5,17],3,29.25]

from txdot_inputs import *
INPUTS = {}
if True:
    for tx in TX:
        #slab thickness
        ts=8.
        ModelName = tx
        y = TX[tx][-1]
        NumberGirders = TX[tx][3]
        GirderSpacing = TX[tx][2]*12.
        SpanLength = TX[tx][1]*12. - 2 * 2
        Overhang = TX[tx][0]*12.
        numberofdiaphragms=TX[tx][5]
        #get girder depth
        Props = {}
        if "R" in tx:
            #rollded girder
            depth = TX[tx][4][0]
            haunch = y-ts-depth
            gz = -ts/2. - haunch - depth/2.
            
            Props[1] = [0, 'AISC','13th','W'+str(TX[tx][4][0])+'x'+str(TX[tx][4][1])]
            Props[2] = [0, 'AISC','13th',DIAPHRAGMS['W'+str(TX[tx][4][0])]]

        elif "P" in tx:
            ttf = TX[tx][4][0]
            btf = TX[tx][4][1]
            tbf = TX[tx][4][2]
            bbf = TX[tx][4][3]
            tw = TX[tx][4][4]
            hw = TX[tx][4][5]
            d = ttf+hw+tbf
            ybar = ComputeSectionProps(ttf,btf,tbf,bbf,tw,hw)
            depth = ttf+tbf+hw
            haunch = y-ts-depth
            gz = -ts/2. - haunch - depth + ybar
            Props[1] = [0, 'IBeam',bbf,btf,d,tbf,ttf,tw]
            Props[2] = [0, 'AISC', '13th',DIAPHRAGMS[d]]

            
        ds = 1./(numberofdiaphragms-1)
        Units = ['in', 'lb', 'Btu', 'psi', 'lbf', 'F', 386.089]
        
        
        SidewalkWidth = 0
        ParapetWidth = 12.
        
        
        DeckThick = ts
        DeckMat = 2
        SidewalkThick = 0.
        ParapetThick = 0.
        ParapetMat = 2
        SidewalkStiffness = False  # thickened shell
        ParapetStiffness = False  # makes a 'girder'
        # LoadLocation = [x,y,fx,fy,fz] #as % of length, width, 0.5,0.5 is center.
        #PointLoad = [[0.5,0.5,0,0,-1000.]]
        PointLoad = []
        #make a few static load cases putting HS20 at midspan, left, right center
        #center
        centerx = SpanLength/2.
        
        totalwidth = (NumberGirders-1)*GirderSpacing+2*Overhang
        shifty = (totalwidth - 2 * ParapetWidth - 24. * 2 - 36. - 2)/2.
        #print totalwidth
        #print shifty
        truckcenters=[(centerx,totalwidth/2.-shifty),(centerx,totalwidth/2.),(centerx,totalwidth/2.+shifty)]
        for tk in truckcenters:
            PointLoad.append([[(tk[0]-224.)/SpanLength,(tk[1]+36)/totalwidth,0,0,-4000.],
                [(tk[0]-224.)/SpanLength,(tk[1]-36)/totalwidth,0,0,-4000.],
                [(tk[0]-56.)/SpanLength, (tk[1]+36)/totalwidth,0,0,-16000.],
                [(tk[0]-56.)/SpanLength, (tk[1]-36)/totalwidth,0,0,-16000.],
                [(tk[0]+112.)/SpanLength,(tk[1]+36)/totalwidth,0,0,-16000.],
                [(tk[0]+112.)/SpanLength,(tk[1]-36)/totalwidth,0,0,-16000.]])
        #-224,-56,112
        #pp.pprint(PointLoad)
        #Vehicles = [['HS20_14FT', 'DESIGN', 0.0]]
        Vehicles=[]
        MPF = [1.2, 1, 0.85, 0.65]
        IM = 1.33
        LF = [1.25,   0.9,     1.5,    0.65,         1.75,         1.35,      1.45,
              1.5,    1.0,    1.0,       1.3,      1.0,         1.3,      1.0]
        PHI = [1, 1]
        # x as % of span length, p as prop integer
        #          [x0,(prop,haunch),x1,(prop,haunch),x2,(prop,haunch),x3]
        Girder_E = [0, (1, haunch), 1.0]
        Girder_I = [0, (1, haunch), 1.0]
        # Girder_common = [MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
        Girder_Icommon = [1, 0, 0, 0, 0, 'C', 1, 1, 0, ds*SpanLength, ds*SpanLength, 0]
        Girder_Ecommon = [1, 0, 0, 0, 0, 'C', 1, 1, 0, ds*SpanLength, ds*SpanLength, 0]
        GirderType = 'C'
        #           [x,mat,prop,z]
        Diaphragm=[]
        for d in range(numberofdiaphragms):
            Diaphragm.append([d*ds, 1, 2, gz])
        Diaphragm_Common = ['straight', 0, 0, 0, 0, 1]



        Mats = {}
        #Mats[no] = [E, nu, rho, strength]
        Mats[1] = [29.0e6, 0.29, 490./12.**3, 50000.]
        Mats[2] = [57000.*(4000.)**0.5, 0.29, 150./12.**3, 4000.]
        deckmesh = 10.
        skew = 0.
        # F for free, R for restrained, float for spring stiffness
        BCS = ['R', 'R', 'R', 'F', 'F', 'F']  # start of each girder
        BCE = ['F', 'F', 'R', 'F', 'F', 'F']  # end of each girder
        BCZ = 0.
        # guidance is truss elements, so release moments
        Options = ['diaphragmreleases', True,
                   'reactions', True, 'displacements', True,'sectioncheckingfigure',False]

        # -------------------------------------------------------------------------------
        # Package them up
        # -------------------------------------------------------------------------------

        variablelist = ['ModelName', 'NumberGirders', 'GirderSpacing', 'SidewalkWidth', 'ParapetWidth', 'SpanLength',
                        'Overhang', 'DeckThick', 'DeckMat', 'SidewalkThick', 'ParapetThick', 'SidewalkStiffness', 'ParapetStiffness',
                        'PointLoad', 'Girder_E', 'Girder_I', 'Diaphragm', 'Props', 'Mats', 'deckmesh', 'skew', 'BCS', 'BCE', 'gz', 'Units',
                        'Options', 'ParapetMat', 'GirderType', 'Girder_Icommon', 'Girder_Ecommon', 'BCZ', 'Diaphragm_Common', 'Vehicles',
                        'MPF', 'IM', 'LF', 'PHI']

        
        INPUTS[tx] = {}
        for var in variablelist:
            INPUTS[tx][var] = eval(var)
        INPUTS[tx]['ModelName'] = tx
        INPUTS[tx]['Label']=tx



# -------------------------------------------------------------------------------
# Parameterize
# -------------------------------------------------------------------------------
if False:
    model = 0
    # SET 1
    # vary the girder spacing
    for i in range(60, 126, 6):
        model += 1
        INPUTS["%03d" % model] = copy.deepcopy(INPUTS['000'])
        INPUTS["%03d" % model]['GirderSpacing'] = i
        INPUTS["%03d" % model]['ModelName'] = ModelName+"%03d" % model
        INPUTS["%03d" % model]['Label'] = 'S='+"%03d" % i+" inches"

if False:
    model = 20
    # SET 2
    # vary the longitudinal bc stiffness
    for k in range(13):
        kstiff = 10**k
        model += 1
        # Start with a wide one, which showed some higher PEM
        INPUTS["%03d" % model] = copy.deepcopy(INPUTS['011'])
        INPUTS["%03d" % model]['BCS'][0] = kstiff
        INPUTS["%03d" % model]['ModelName'] = ModelName+"%03d" % model
        INPUTS["%03d" % model]['Label'] = 'BC_kx='+"%4.1e" % kstiff+" lbf/in"

if False:
    # SET 3
    # vary the span length? (set 1 showed higher PEM for larger beam spacing, but i suspect that is actually the result of bridge plan aspect ratio, test that in set 3)
    model = 40
    # vary the longitudinal bc stiffness
    for k in range(11):
        model += 1
        # Start with a wide one, which showed some higher PEM
        INPUTS["%03d" % model] = copy.deepcopy(INPUTS['011'])
        INPUTS["%03d" % model]['SpanLength'] = 780.*(0.5+(k/10.))
        INPUTS["%03d" % model]['ModelName'] = ModelName+"%03d" % model
        INPUTS["%03d" % model]['Label'] = 'Span=' + \
            "%03d" % (INPUTS["%03d" % model]['SpanLength']/12.)+" ft"

if False:
    # SET 4, one without deck shear stiffness
    # this was BCW57901_033 & _034 which showed the big difference in deflection as a result of deck shear modulus
    # SET 2 didn't show much difference in deflection so maybe deflection is a result of G and not PEM.
    model = 60
    if True:
        model += 1
        # Start with a wide one, which showed some higher PEM
        INPUTS["%03d" % model] = copy.deepcopy(INPUTS['011'])
        INPUTS["%03d" % model]['Mats'][2] = [57000.*(4000.)**0.5, 57000.*(4000.)**0.5, 57000.*(4000.)**0.5, 0.0, 0.0, 0.0, 150./12.**3, 4000., (57000.*(
            4000.)**0.5 / (2*(1+0.2)))*0.000001, 57000.*(4000.)**0.5 / (2*(1+0.2)), 57000.*(4000.)**0.5 / (2*(1+0.2))]
        INPUTS["%03d" % model]['ModelName'] = ModelName+"%03d" % model
        INPUTS["%03d" % model]['Label'] = 'Zero Deck Shear Modulus'

#print len(INPUTS)

# -------------------------------------------------------------------------------
# Make the Runs
# -------------------------------------------------------------------------------
#run1 = [inp for inp in sorted(INPUTS.keys()) if int(inp) >= 0 and int(inp) < 20]
#run2 = [inp for inp in sorted(INPUTS.keys()) if int(inp) > 20 and int(inp) < 40]
#run3 = [inp for inp in sorted(INPUTS.keys()) if int(inp) > 40 and int(inp) < 60]
#run4 = [inp for inp in sorted(INPUTS.keys()) if int(inp) > 60 and int(inp) < 80] + ['011']
run1 = sorted(INPUTS.keys())
run = copy.deepcopy(run1) 

if True:
    for inp in run1:
        WriteInput(INPUTS[inp])

if False:
    for inp in run2:
        WriteInput(INPUTS[inp])

    for inp in run3:
        WriteInput(INPUTS[inp])

    for inp in run4:
        WriteInput(INPUTS[inp])
run = copy.deepcopy(run1) 

#CHANGE THIS TO USE MULTI-PROESSING, 6 AT A TIME?
if True:
    # first is use to run 423, after that, runs are submitted with RR3.
    First = True

    # NEWRUNS = True will force new runs for each inp, NEWRUNS=False will attempt to use existing, and only rerun if necessary.
    NEWRUNS = False
    while len(run) > 0:
        rerun = []
        for inp in run:
            success = False
            print "Running "+inp
            if NEWRUNS:
                if First:
                    # WriteInput(INPUTS[inp])
                    #RunMiser(INPUTS[inp], 4, 2, 3)
                    RunMiser(INPUTS[inp], 3, 1, 1)
                else:
                    RunMiser(INPUTS[inp], 'R', 'R', 3)

            # check that it completed successfully
            logfile = GetLogFile(INPUTS[inp]['ModelName'])
            if os.path.isfile(logfile):
                with open(logfile, 'r') as f:
                    if "ELAPSED SOLUTION TIME" in f.readlines()[-2]:
                        success = True
                        print "Success"
                    if not success:
                        print "Failed"
                        # put it back in the queue
                        if NEWRUNS:
                            rerun.append(inp)
                        else:
                            RunMiser(INPUTS[inp], 'R', 'R', 1)
            else:
                # WriteInput(INPUTS[inp])
                #RunMiser(INPUTS[inp], 4, 2, 3)
                RunMiser(INPUTS[inp], 3, 1, 1)

        if First:
            First = False
        run = copy.deepcopy(rerun)

    run = copy.deepcopy(run1) 

# -------------------------------------------------------------------------------
# Fetch the model data
# -------------------------------------------------------------------------------
# open the shelves and grab:
    # the elements & nodes along each girder line
    # the node coords of BCs
MODELDATA = {}
recoverkeys = ['GIRDERLINEELEMENTS', 'NODES', 'BEAMS', 'BCnodes', 'G_listofelements',
               'NODES_OUT', 'G_listofys', 'COMP', 'NONCOMP', 'REACTIONS', 'DISPLACEMENTS', 'NODEMAP']

if False:
    for inp in sorted(run1):
        MODELDATA[inp] = recovershelve(
            INPUTS[inp]['ModelName'], [], recoverkeys)
    for inp in sorted(run2):
        MODELDATA[inp] = recovershelve(
            INPUTS[inp]['ModelName'], [], recoverkeys)
    for inp in sorted(run3):
        MODELDATA[inp] = recovershelve(
            INPUTS[inp]['ModelName'], [], recoverkeys)
    # for inp in sorted(run):
    #    MODELDATA[inp] = recovershelve(INPUTS[inp]['ModelName'],[],recoverkeys)
    for inp in sorted(run4):
        MODELDATA[inp] = recovershelve(
            INPUTS[inp]['ModelName'], [], recoverkeys)
if True:
    for inp in sorted(run):
        MODELDATA[inp] = recovershelve(
            INPUTS[inp]['ModelName'], [], recoverkeys)

# -------------------------------------------------------------------------------
# Plot Results
# -------------------------------------------------------------------------------
if False:

    MomentPlot(MODELDATA, INPUTS, run1,
               'ParametricController_MomentSummary1', norm=True)
    MomentPlot(MODELDATA, INPUTS, run2,
               'ParametricController_MomentSummary2', norm=True)
    MomentPlot(MODELDATA, INPUTS, run3,
               'ParametricController_MomentSummary3', norm=True)

    DeflectionPlot(MODELDATA, INPUTS, run1,
                   'ParametricController_DeflectionSummary1')
    DeflectionPlot(MODELDATA, INPUTS, run2,
                   'ParametricController_DeflectionSummary2')
    DeflectionPlot(MODELDATA, INPUTS, run3,
                   'ParametricController_DeflectionSummary3')

    MomentPlot(MODELDATA, INPUTS, run4,
               'ParametricController_MomentSummary4', norm=True)
    DeflectionPlot(MODELDATA, INPUTS, run4,
                   'ParametricController_DeflectionSummary4')
if False:

    MomentPlot(MODELDATA, INPUTS, run,
               'ParametricController_MomentSummary_TX', norm=True)
    DeflectionPlot(MODELDATA, INPUTS, run,
                   'ParametricController_DeflectionSummary_TX')

#build dataframe of results
if True:
    RESULTS= {}
    for inp in run:
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
                minendmoment = min(minendmoment, (float(max([MODELDATA[inp]['COMP'][g][section][c][4] for c in cases]))))

        RESULTS[inp]['maxendmoment'] = maxendmoment
        RESULTS[inp]['minendmoment'] = minendmoment

        

resDF = pd.DataFrame.from_dict(RESULTS,orient='index')
resDF.to_csv('txdot_results.csv')