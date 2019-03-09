


from subprocess import CREATE_NEW_CONSOLE
import subprocess
from subprocess import STDOUT
import copy
import os
#import pandas as pd
import math
import re
import shelve
import plotly.offline as py
import plotly.graph_objs as go

def ScrubLogFile(logfilename,masssum=False):
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
    LogRead=open(logfilename,'r')
    LineRead = LogRead.readline()
    solution = 'fail'
    numberofnotes='X'
    numberofwarnings='X'
    numberoferrors='X'
    FZList=[]
    TimeList=[]
    try:
        while LineRead:
            if LineRead.split(":")[0].strip() == "*Number of Notes":
                numberofnotes=LineRead.split(":")[1].strip()
            if LineRead.split(":")[0].strip() == "*Number of Warnings":
                numberofwarnings=LineRead.split(":")[1].strip()
            if LineRead.split(":")[0].strip() == "*Number of Errors":
                numberoferrors=LineRead.split(":")[1].strip()
            if "*Solution terminated" in LineRead:
                solution='Terminated'        
            if "*Solution completed" in LineRead:
                solution='Completed'
            if "DIRECT SUMMATION OF NODE REACTION FORCES" in LineRead:
                #print LineRead.split()
                try:
                    timestep = int(LineRead.split()[8])
                except IndexError:
                    timestep = 0
                LineRead = LogRead.readline()
                LineRead = LogRead.readline()
                FZ=float(LineRead.split()[3])
                FZList.append(FZ)
                TimeList.append(timestep)
            LineRead = LogRead.readline()
        LogRead.close()
    except:
        pass
    #print logfilename,TimesOfLessThanMaxReaction
    if masssum:
        return numberofnotes,numberofwarnings,numberoferrors,solution,maxreaction
    else:
        return numberofnotes,numberofwarnings,numberoferrors,solution


def GetLogFile(OutputFile):
    path = os.getcwd()
    files = os.listdir(path)
    logfiles = []
    #add date regex
    pattern1=re.compile("[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]")

    for f in files:
        if f[-3:]=='log' and OutputFile in f and bool(re.search(pattern1,f)) and 'logger' not in f:
            logfiles.append(f)
    if logfiles == []:
        if os.path.isfile(OutputFile+".shv"):
            logfiles.append(OutputFile+".shv")
    if len(logfiles)==0:
        return 'None'
    else:
        logfile = sorted(logfiles)[-1]
        return logfile


def recovershelve(OutputFile,skips=[],only=False):
    '''
    function to recover a shelve contents
    outputfile is the name of the shelf (without date & time)
    skips are the contents to skip 
    only can be specified to collect only certain variables
    '''
    #find the most recent shelve file...
    path = os.getcwd()
    files = os.listdir(path)
    shelves=[]

    #add date regex
    pattern1=re.compile("[0-9][0-9][0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]_[0-9][0-9]")

    for f in files:
        if f[-3:]=='shv' and OutputFile in f and bool(re.search(pattern1,f)):
            shelves.append(f)
    if shelves == []:
        if os.path.isfile(OutputFile+".shv"):
            shelves.append(OutputFile+".shv")

    shelffile = sorted(shelves)[-1]
    #print "recovering "+shelffile
    RETURNS={}

    my_shelf = shelve.open(shelffile,'r')
    for key in sorted(my_shelf):
        if not only:
            if key not in ['bridgemiserstart','timepoint1','timepoint0','timeelapsed','Logger','Tee','color_cycle','FILE','overwriteshelveratinginputs','SKIPSHELVES','DEBUG_MODE','DEBUGFILE','Loggertimestring'] and key not in skips:
                try:
                    RETURNS[key]=my_shelf[key]
                except:
                    printred("Cannot unshelve key: "+key)
        else:
            if key in only:
                try:
                    RETURNS[key]=my_shelf[key]
                    #print key
                except:
                    printred("Cannot unshelve key: "+key)
    my_shelf.close()
    return RETURNS


#-------------------------------------------------------------------------------
# Write input file
#-------------------------------------------------------------------------------
def WriteInput(INPUTS):
    #now build the input file(s)
    inpfilename = INPUTS['ModelName']+".inp"
    inpfile = open(inpfilename,'w')

    inpfile.write("*OUTPUTFILE,"+INPUTS['ModelName']+",\n")
    inpfile.write("*UNITS,"+",".join([str(j) for j in INPUTS['Units']])+",\n")
    inpfile.write("*OPTIONS,"+",".join([str(j) for j in INPUTS['Options']])+",\n")
    inpfile.write("*DECKMESH,"+str(INPUTS['deckmesh'])+",\n")
    for m in sorted(INPUTS['Mats'].keys()):
        inpfile.write("*MATERIAL,"+str(int(m))+","+",".join([str(j) for j in INPUTS['Mats'][m]])+",\n")
    for bp in sorted(INPUTS['Props'].keys()):
        inpfile.write("*BEAMPROP,"+str(int(bp))+","+",".join([str(j) for j in INPUTS['Props'][bp]])+",\n")

    #DECK-------------------------------------------------------------------
    d1W = (INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] - 1) + 2 * INPUTS['Overhang']) - 2 * INPUTS['SidewalkWidth']
    d1L = INPUTS['SpanLength']

    #assumed D1 is at z=0
    DL0=0
    DL1=0
    DL2=0
    DW=0
    KF=0
    inpfile.write("*DECK,1,"+str(d1L)+','+str(d1W)+',0,'+str(INPUTS['DeckMat'])+','+str(INPUTS['DeckThick'])+',0,'+str(-d1W/2.)+','+','.join([str(j) for j in [DL0,DL1,DL2,DW,KF]])+"\n")

    d2W = INPUTS['SidewalkWidth']
    d2L = INPUTS['SpanLength']

    if INPUTS['SidewalkThick']>0 and d2W>0 and INPUTS['SidewalkStiffness']:
        #ADD DECK 2 and 3 with thicker deck
        z = (INPUTS['DeckThick'] + INPUTS['SidewalkThick']) / 2. - INPUTS['DeckThick'] / 2.
        inpfile.write("*DECK,2,"+str(d2L)+','+str(d2W)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(INPUTS['DeckThick']+INPUTS['SidewalkThick'])+',0,'+str(-d1W/2. - d2W)+','+','.join([str(j) for j in [DL0,DL1,DL2,DW,KF]])+"\n")
        inpfile.write("*DECK,3,"+str(d2L)+','+str(d2W)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(INPUTS['DeckThick']+INPUTS['SidewalkThick'])+',0,'+str(d1W/2.)       +','+','.join([str(j) for j in [DL0,DL1,DL2,DW,KF]])+"\n")

    elif not INPUTS['SidewalkStiffness']:
        #ADD DECK 2 and 3 with NSM in SIDL
        z = 0.
        DW = INPUTS['SidewalkThick'] * INPUTS['Mats'][INPUTS['DeckMat']][2]
        inpfile.write("*DECK,2,"+str(d2L)+','+str(d2W)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(INPUTS['DeckThick'])+',0,'+str(-d1W/2. - d2W)+','+','.join([str(j) for j in [DL0,DL1,DL2,DW,KF]])+"\n")
        inpfile.write("*DECK,3,"+str(d2L)+','+str(d2W)+','+str(z)+','+str(INPUTS['DeckMat'])+','+str(INPUTS['DeckThick'])+',0,'+str(d1W/2.)       +','+','.join([str(j) for j in [DL0,DL1,DL2,DW,KF]])+"\n")

    #GIRDERS-------------------------------------------------------------------
    gx=0
    y0 = -(INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] - 1))/2.
    BCs=[]
    BCe=[]
    for g in range(INPUTS['NumberGirders']):
        #*GIRDER,No,X1,X2,Y,Haunch,Prop,Mat,DL0,DL1,DL2,DW,Comp,Deck,KF,Tfubl,Bfubl,TSS,Cope
        if g in [0,INPUTS['NumberGirders']-1]:
            #exterior
            for x in range(int((len(INPUTS['Girder_E'])-1)/2.)):
                gx+=1        
                xstart = INPUTS['Girder_E'][2*x]*d1L
                prop = INPUTS['Girder_E'][2*x+1][0]
                haunch = INPUTS['Girder_E'][2*x+1][1]
                xend = INPUTS['Girder_E'][2*x+2]*d1L
                inpfile.write("*GIRDER,"+str(gx)+','+str(xstart)+","+str(xend)+','+str(y0)+','+str(haunch)+','+str(prop))
                #Girder_common = [haunch, MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
                inpfile.write(','+','.join([str(j) for j in INPUTS['Girder_Ecommon']])+"\n")
                if x==0:
                    BCs.append(gx)
                elif x==int((len(INPUTS['Girder_E'])-1)/2.)-1:
                    BCe.append(gx)
        else:
            #interior
            for x in range(int((len(INPUTS['Girder_I'])-1)/2.)):
                gx+=1        
                xstart = INPUTS['Girder_I'][2*x]*d1L
                prop = INPUTS['Girder_I'][2*x+1][0]
                haunch = INPUTS['Girder_I'][2*x+1][1]
                xend = INPUTS['Girder_I'][2*x+2]*d1L
                inpfile.write("*GIRDER,"+str(gx)+','+str(xstart)+","+str(xend)+','+str(y0)+','+str(haunch)+','+str(prop))
                #Girder_common = [haunch, MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
                inpfile.write(','+','.join([str(j) for j in INPUTS['Girder_Icommon']])+"\n")
                if x==0:
                    BCs.append(gx)
                elif x==int((len(INPUTS['Girder_E'])-1)/2.)-1:
                    BCe.append(gx)
        y0+=INPUTS['GirderSpacing']

    #compute deck inputs
    if INPUTS['ParapetWidth']>0 and INPUTS['ParapetThick']>0:
        #make a new beam prop
        inpfile.write("*BEAMPROP,"+str(int(999))+',0,SolidTrapezoid,'+str(INPUTS['ParapetWidth'])+','+str(INPUTS['ParapetWidth'])+','+str(INPUTS['ParapetThick'])+'\n')
        #now the girders
        yp = ((INPUTS['GirderSpacing'] * (INPUTS['NumberGirders'] - 1) + 2 * INPUTS['Overhang']) - INPUTS['ParapetWidth'])/2.
        haunch = -INPUTS['DeckThick']-INPUTS['ParapetThick']
        if INPUTS['SidewalkThick']>0 and INPUTS['SidewalkWidth']>0 and INPUTS['SidewalkStiffness']:
            haunch -= INPUTS['SidewalkThick']
        inpfile.write("*GIRDER,998,0,"+str(d1L)+','+str(-yp)+','+str(haunch)+',999,'+str(INPUTS['DeckMat'])+',0,0,0,0,B,1,3\n')
        inpfile.write("*GIRDER,999,0,"+str(d1L)+','+str(yp)+','+str(haunch)+',999,'+str(INPUTS['DeckMat']) +',0,0,0,0,B,1,3\n')
        #*GIRDER, Girder Number, Xstart, Xend, Ycg, haunch, BeamPropNo, MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP (C;NC;NA), corresponding deck no, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE

    
    #*BC_GIRDER, Girder Number, X, Z, DOF1, DOF2, DOF3, DOF4, DOF5, DOF6 ,,,,,,,,,
    for gs in BCs:
        inpfile.write("*BC_GIRDER,"+str(gs)+',START,'+str(INPUTS['BCZ'])+','+','.join([str(j) for j in INPUTS['BCS']])+"\n")
    for ge in BCe:
        inpfile.write("*BC_GIRDER,"+str(ge)+',END,'+str(INPUTS['BCZ'])+','+','.join([str(j) for j in INPUTS['BCE']])+"\n")
        

    #DIAPHRAGMS-------------------------------------------------------------------
    #*DIAPHRAGM,No, zi, zm, zj,Mat,Prop,X, type,DL0,DL1,DL2,DW,KF,skip,,,,
    dno=0
    for d in INPUTS['Diaphragm']:
        dno+=1
        inpfile.write('*DIAPHRAGM,'+str(dno)+','+str(d[3])+','+str(d[3])+','+str(d[3])+','+str(d[1])+','+str(d[2])+','+str(d[0]*d1L)+',')
        inpfile.write(','.join([str(j) for j in INPUTS['Diaphragm_Common']])+'\n')
    
    if len(INPUTS['PointLoad'])>0:
        inpfile.write('*STATICLOADCASE,100,PointLoad\n')
        for pt in INPUTS['PointLoad']:
            inpfile.write('*NODEFORCE,'+str(pt[0]*d1L)+','+str(pt[1]*d1W - d1W/2.)+',0,100,'+str(pt[2])+','+str(pt[3])+','+str(pt[4])+',0,0,0\n')
    
    if len(INPUTS['Vehicles'])>0:
        #*LANE, Lane Number, shape, curvature, x1, x2, y, z, width
        #Fit lanes between sidewalks
        if d1W<=17*12:
            numlanes = 1
            if d1W<=144:
                #one lane
                width = d1W
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L)+',0,0,'+str(width-20.)+'\n')
            else:
                width = 144
                #-
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
                #0
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L)+',0,0,'+str(width)+'\n')
                #+
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANESETS, 1\n')
                inpfile.write('*LANESETS, 2\n')
                inpfile.write('*LANESETS, 3\n')

        elif d1W<=308:
            numlanes = 2
            width = (d1W-48.)/2.

            #0
            inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
            inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')
        else:
            numlanes = math.floor(d1W/144.)
            width = 144.-24.
            if numlanes == 2:
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+3*width/2.+48.)+',0,'+str(width)+'\n')
                #0
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-width/2.-24)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,4,STRAIGHT,FLAT,0,'+str(d1L)+','+str(width/2.+24)+',0,'+str(width)+'\n')
                #+
                inpfile.write('*LANE,5,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-3*width/2.-48.)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,6,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANESETS, 1,2\n')
                inpfile.write('*LANESETS, 3,4\n')
                inpfile.write('*LANESETS, 5,6\n')

            elif numlanes == 3:
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+width/2.+5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+3*width/2.+48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+5*width/2.+96)+',0,'+str(width)+'\n')
                #0
                inpfile.write('*LANE,4,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-2*width/2.-24)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,5,STRAIGHT,FLAT,0,'+str(d1L)+','+str(0)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,6,STRAIGHT,FLAT,0,'+str(d1L)+','+str(2*width/2.+24)+',0,'+str(width)+'\n')
                #+
                inpfile.write('*LANE,7,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-5*width/2.-96)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,8,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-3*width/2.-48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,9,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-width/2.-5)+',0,'+str(width-10)+'\n')

                inpfile.write('*LANESETS, 1,2,3\n')
                inpfile.write('*LANESETS, 4,5,6\n')
                inpfile.write('*LANESETS, 7,8,9\n')                
            else:
                if numlanes > 4:
                    print "limiting lanes to 4"
                    numlanes == 4
                inpfile.write('*LANE,1,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+geomwidth/2.+5)+',0,'+str(width-10)+'\n')
                inpfile.write('*LANE,2,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+3*geomwidth/2.+48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,3,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+5*geomwidth/2.+96)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,4,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-d1W/2.+7*geomwidth/2.+144)+',0,'+str(width)+'\n')
                #0
                inpfile.write('*LANE,5,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-3*geomwidth/2.-72)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,6,STRAIGHT,FLAT,0,'+str(d1L)+','+str(-geomwidth/2.-24)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,7,STRAIGHT,FLAT,0,'+str(d1L)+','+str(geomwidth/2.+24)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,8,STRAIGHT,FLAT,0,'+str(d1L)+','+str(3*geomwidth/2.+72)+',0,'+str(width)+'\n')
                #+
                inpfile.write('*LANE,9,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-7*geomwidth/2.-144)+',0,'+str(width)+'\n')
                inpfile.write('*LANE,10,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-5*geomwidth/2.-96)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,11,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-3*geomwidth/2.-48)+',0,'+str(width-24)+'\n')
                inpfile.write('*LANE,12,STRAIGHT,FLAT,0,'+str(d1L)+','+str(d1W/2.-geomwidth/2.-5)+',0,'+str(width-10)+'\n')

                inpfile.write('*LANESETS, 1,2,3,4\n')
                inpfile.write('*LANESETS, 5,6,7,8\n')
                inpfile.write('*LANESETS, 9,10,11,12\n')


        for v in INPUTS['Vehicles']:
            inpfile.write('*VEHICLE,20,'+str(v[0])+','+str(v[1])+','+str(v[2])+'\n')
        inpfile.write('*MPF,'+','.join([str(j) for j in INPUTS['MPF']])+'\n')
        inpfile.write('*IMPACTFACTOR,'+str(INPUTS['IM'])+'\n')
        inpfile.write('*LOADFACTOR,'+','.join([str(j) for j in INPUTS['LF']])+'\n')
        inpfile.write('*PHI,'+','.join([str(j) for j in INPUTS['PHI']])+'\n')

    inpfile.write('*PLOTLYMESH\n')
    inpfile.write('*PLOTLYSOLID\n')
    inpfile.write('*MULTI,6,6\n')


    inpfile.close()


#-------------------------------------------------------------------------------
# Run Miser
#-------------------------------------------------------------------------------
def RunMiser(INPUTS,runMODE,extractMODE,rateMODE,skipshelves=False):

    FileName = INPUTS['ModelName']+'.inp'
    try:
        if skipshelves:
            a = subprocess.check_call(["python32",'bridgemiser_016.py',FileName,str(runMODE),str(extractMODE),str(rateMODE),'on','on'],stderr=STDOUT)
        else:
            a = subprocess.check_call(["python32",'bridgemiser_016.py',FileName,str(runMODE),str(extractMODE),str(rateMODE)],stderr=STDOUT)
    except Exception as e:
        #raise e
        return False
    return True



def MomentPlot(MODELDATA,INPUTS,run,filename,norm=False):
    tableau10_256=[]
    tableau20=[(31,119,180),(174,199,232),(255,127,14),(255,187,120),(44,160,44),(152,223,138),(214,39,40),(255,152,150),(148,103,189),(197,176,213),(140,86,75),(196,156,148),(227,119,194),(247,182,210),(127,127,127),(199,199,199),(188,189,34),(219,219,141),(23,190,207),(158,218,229)]
    for i in range(len(tableau20)/2):   #integer operation
        r, g, b = tableau20[i*2]
        #color_cycle2.append((r / 255., g / 255., b / 255.))
        tableau10_256.append((r, g, b))


    PlotTraces=[]
    maxz=0
    maxy=0
    miny=0
    maxx=0
    minx=0
    for inp in sorted(run):
        first=True
        cc = sorted(INPUTS.keys()).index(inp)%10
        maxZinp=0
        if norm:
            for yg in sorted(MODELDATA[inp]['GIRDERLINEELEMENTS'].keys()):
                for tup in MODELDATA[inp]['GIRDERLINEELEMENTS'][yg]:
                    g = 'g'+str(tup[2])
                    section="SC_"+"%04d"%tup[0]
                    maxZinp = max(maxZinp,(float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4])))

        for yg in sorted(MODELDATA[inp]['GIRDERLINEELEMENTS'].keys()):
            x=[]
            y=[]
            zmin=[]
            zmax=[]
            hovermax=[]
            hovermin=[]
            for tup in MODELDATA[inp]['GIRDERLINEELEMENTS'][yg]:
                #list of (bee,xi,g,xj)
                x.append(float(tup[1]))
                y.append(float(yg))
                g = 'g'+str(tup[2])
                section="SC_"+"%04d"%tup[0]
                x.append(float(tup[3]))
                y.append(float(yg))
                if norm:
                    zmax.append(float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4])/maxZinp)
                    zmin.append(float(MODELDATA[inp]['COMP'][g][section][6]['MIN'][4])/maxZinp)
                    zmax.append(float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4+6])/maxZinp)
                    zmin.append(float(MODELDATA[inp]['COMP'][g][section][6]['MIN'][4+6])/maxZinp)
                    hovermax.append(float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4]))
                    hovermin.append(float(MODELDATA[inp]['COMP'][g][section][6]['MIN'][4]))
                    hovermax.append(float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4+6]))
                    hovermin.append(float(MODELDATA[inp]['COMP'][g][section][6]['MIN'][4+6]))
                else:
                    zmax.append(float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4])/maxZinp)
                    zmin.append(float(MODELDATA[inp]['COMP'][g][section][6]['MIN'][4])/maxZinp)
                    zmax.append(float(MODELDATA[inp]['COMP'][g][section][6]['MAX'][4+6])/maxZinp)
                    zmin.append(float(MODELDATA[inp]['COMP'][g][section][6]['MIN'][4+6])/maxZinp)
                #compoutputterms = ['Axial','HShear','VShear','Torsion','HMoment','VMoment']
            maxz = max(maxz,max(zmax))
            maxy = max(maxy,max(y))
            miny = min(miny,min(y))
            maxx = max(maxx,max(x))
            minx = min(minx,min(x))
            if norm:
                if first:
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmax,hoverinfo="all",text=hovermax,legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly'))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmin,hoverinfo="all",text=hovermin,legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly',showlegend=False))
                    first=False
                else:
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmax,hoverinfo="all",text=hovermax,legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly',showlegend=False))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmin,hoverinfo="all",text=hovermin,legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly',showlegend=False))
            else:
                if first:
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmax,hoverinfo="all",legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly'))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmin,hoverinfo="all",legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly',showlegend=False))
                    first=False
                else:
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmax,hoverinfo="all",legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly',showlegend=False))
                    PlotTraces.append(dict(type='scatter3d', line=dict(color='rgb'+`tableau10_256[cc]`,width=3), mode='lines',x=x,y=y,z=zmin,hoverinfo="all",legendgroup=inp,name=INPUTS[inp]['ModelName']+" "+INPUTS[inp]['Label'],visible='legendonly',showlegend=False))

    noaxis=dict(showbackground=False,showline=False,zeroline=False,showgrid=False,showticklabels=False,title='')
    layout = dict(
        title="Moment Summary",showlegend=True,
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(x=(maxx-minx)/(maxy-miny), y=1, z=1),
            camera=dict(eye=dict(x=0, y=-0.5, z=0)),
            xaxis=noaxis,
            yaxis=noaxis,
            zaxis=noaxis))
    fig = dict(data=PlotTraces,layout=layout)
    py.plot(fig, filename=filename+'.html',auto_open=False,show_link=False)


#-------------------------------------------------------------------------------
# Import Results
#-------------------------------------------------------------------------------
#def ReadResults(INPUTS):
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

#-------------------------------------------------------------------------------
# INPUTS
#-------------------------------------------------------------------------------
if True:
    ModelName = 'MBEA1_P'
    gz=-24.0
    Units = ['in','lb','Btu','psi','lbf','F',386.089]
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
    SidewalkStiffness = True   #thickened shell
    ParapetStiffness = True    #makes a 'girder'
    #LoadLocation = [x,y,fx,fy,fz] #as % of length, width, 0.5,0.5 is center.
    #PointLoad = [[0.5,0.5,0,0,-1000.]]
    PointLoad = []
    Vehicles=[['HS20_14FT','DESIGN',0.0]]
    MPF = [1.2,1,0.85,0.65]
    IM = 1.33
    LF = [1.25,   0.9,     1.5,    0.65,         1.75,         1.35,      1.45,      1.5,    1.0,    1.0,       1.3,      1.0,         1.3,      1.0]
    PHI=[1,1]
    #x as % of span length, p as prop integer
    #          [x0,(prop,haunch),x1,(prop,haunch),x2,(prop,haunch),x3]
    Girder_E = [0,(1,0),0.2,(2,0),0.8,(1,0),1.0]
    Girder_I = [0,(1,0),0.2,(3,0),0.8,(1,0),1.0]
    #Girder_common = [MaterialNo, DL0-NSM, DL1-NSM, DL2-NSM, SIDL-NSM, COMP, Stiffness Factor, Top Flange Unbraced Length, Bottom Flange Unbraced Length, transverse stiffener spacing, COPE]
    Girder_Icommon=[1,0,0,0,0,'C',1,1,0,195,195,0]
    Girder_Ecommon=[1,0,0,0,0,'C',1,1,0,195,195,0]
    GirderType = 'C'
    #           [x,mat,prop,z]
    Diaphragm = [[0,1,4,gz],[0.25,1,4,gz],[0.5,1,4,gz],[0.75,1,4,gz],[1,1,4,gz]]
    Diaphragm_Common=['straight',0,0,0,0,1]
    Props={}
    Props[1]=[0,'IBeam',11.5,11.5,33.1,0.855,0.855,0.58]
    Props[2]=[0,'IBeam',11.0835,11.5,33.85,1.48,0.855,0.58]
    Props[3]=[0,'IBeam',11.038,11.5,33.725,1.605,0.855,0.58]
    Props[4]=[0,'Lipchannel',3.95,18,0,0.625,0.45,0]

    #Props[1]=[1,0,Null,37.8712,6607.75,217.235,6.73246,40.75,20.375]
    #Props[2]=[2,0,Null,44.5148,8277.46,276.798,15.5185,40.75,20.375]
    #Props[3]=[3,0,Null,45.6822,8414.14,288.743,18.3722,40.75,20.375]
    #Props[4]=[4,0,Null,12.475,549.035,15.6833,1.10832,18,9]
    Mats={}
    #Mats[no] = [E, nu, rho, strength]
    Mats[1] = [29.0e6, 0.29, 490./12.**3, 36000.]
    Mats[2] = [57000.*(4000.)**0.5, 0.29, 150./12.**3, 4000.]
    deckmesh = 10.
    skew = 0.
    #F for free, R for restrained, float for spring stiffness
    BCS = ['R','R','R','F','F','F'] #start of each girder
    BCE = ['F','F','R','F','F','F'] #end of each girder
    BCZ = 0.
    #guidance is truss elements, so release moments
    Options=['diaphragmreleases',True,'reactions',True,'displacements',True]



    #-------------------------------------------------------------------------------
    # Package them up
    #-------------------------------------------------------------------------------

    variablelist = ['ModelName','NumberGirders','GirderSpacing','SidewalkWidth','ParapetWidth','SpanLength',
    'Overhang','DeckThick','DeckMat','SidewalkThick','ParapetThick','SidewalkStiffness','ParapetStiffness',
    'PointLoad','Girder_E','Girder_I','Diaphragm','Props','Mats','deckmesh','skew','BCS','BCE','gz','Units',
    'Options','ParapetMat','GirderType','Girder_Icommon','Girder_Ecommon','BCZ','Diaphragm_Common','Vehicles',
    'MPF','IM','LF','PHI']

    INPUTS={}
    INPUTS['000']={}
    for var in variablelist:
        INPUTS['000'][var]=eval(var)
    INPUTS['000']['ModelName'] = ModelName+'000'


#-------------------------------------------------------------------------------
# Parameterize
#-------------------------------------------------------------------------------
model=0
#SET 1 
# vary the girder spacing
for i in range(60,126,6):
    model+=1
    INPUTS["%03d"%model] = copy.deepcopy(INPUTS['000'])
    INPUTS["%03d"%model]['GirderSpacing']=i
    INPUTS["%03d"%model]['ModelName'] = ModelName+"%03d"%model
    INPUTS["%03d"%model]['Label'] = 'S='+"%03d"%i+" inches"



model=20
#SET 2
# vary the longitudinal bc stiffness
for k in range(13):
    kstiff = 10**k
    model+=1
    #Start with a wide one, which showed some higher PEM
    INPUTS["%03d"%model] = copy.deepcopy(INPUTS['011'])
    INPUTS["%03d"%model]['BCS'][0]=kstiff
    INPUTS["%03d"%model]['ModelName'] = ModelName+"%03d"%model
    INPUTS["%03d"%model]['Label'] = 'BC_kx='+"%4.1e"%kstiff+" lbf/in"


#SET 3
# vary the span length? (set 1 showed higher PEM for larger beam spacing, but i suspect that is actually the result of bridge plan aspect ratio, test that in set 3)

model=40
# vary the longitudinal bc stiffness
for k in range(11):
    model+=1
    #Start with a wide one, which showed some higher PEM
    INPUTS["%03d"%model] = copy.deepcopy(INPUTS['011'])
    INPUTS["%03d"%model]['SpanLength']= 780.*(0.5+(k/10.))
    INPUTS["%03d"%model]['ModelName'] = ModelName+"%03d"%model
    INPUTS["%03d"%model]['Label'] = 'Span='+"%03d"%(INPUTS["%03d"%model]['SpanLength']/12.)+" ft"








#-------------------------------------------------------------------------------
# Make the Runs 
#-------------------------------------------------------------------------------
run1 = [inp for inp in sorted(INPUTS.keys()) if int(inp)>0 and int(inp)<20]
run2 = [inp for inp in sorted(INPUTS.keys()) if int(inp)>20 and int(inp)<40]
run3 = [inp for inp in sorted(INPUTS.keys()) if int(inp)>40]


run=copy.deepcopy(run3)

if False:
    for inp in run1:
        WriteInput(INPUTS[inp])

    for inp in run2:
        WriteInput(INPUTS[inp])

    for inp in run3:
        WriteInput(INPUTS[inp])


if True:
    #first is use to run 423, after that, runs are submitted with RR3.
    First=True

    #NEWRUNS = True will force new runs for each inp, NEWRUNS=False will attempt to use existing, and only rerun if necessary. 
    NEWRUNS =  False
    while len(run)>0:
        rerun=[]
        for inp in run:
            success = False
            print inp
            if NEWRUNS:
                if First:
                    #WriteInput(INPUTS[inp])
                    RunMiser(INPUTS[inp],4,2,3)
                else:
                    RunMiser(INPUTS[inp],'R','R',3)

            #check that it completed successfully
            logfile = GetLogFile(INPUTS[inp]['ModelName'])
            if os.path.isfile(logfile):
                with open(logfile,'r') as f:
                    if "ELAPSED SOLUTION TIME" in f.readlines()[-2]:
                        success = True
                        print "Success"
                    if not success:
                        print "Failed"
                        #put it back in the queue
                        if NEWRUNS:
                            rerun.append(inp)
                        else:
                            RunMiser(INPUTS[inp],'R','R',3)
            else:
                #WriteInput(INPUTS[inp])
                RunMiser(INPUTS[inp],4,2,3)

        if First:
            First=False
        run = copy.deepcopy(rerun)

#-------------------------------------------------------------------------------
# Fetch the model data
#-------------------------------------------------------------------------------
if True:
    #open the shelves and grab:
        # the elements & nodes along each girder line
        # the node coords of BCs
    MODELDATA={}
    recoverkeys = ['GIRDERLINEELEMENTS','NODES','BEAMS','BCnodes','G_listofelements','NODES_OUT','G_listofys','COMP','NONCOMP','REACTIONS','DISPLACEMENTS']
    #for inp in sorted(run1):
    #    MODELDATA[inp] = recovershelve(INPUTS[inp]['ModelName'],[],recoverkeys)
    #for inp in sorted(run2):
    #    MODELDATA[inp] = recovershelve(INPUTS[inp]['ModelName'],[],recoverkeys)
    #for inp in sorted(run3):
    #    MODELDATA[inp] = recovershelve(INPUTS[inp]['ModelName'],[],recoverkeys)
    for inp in sorted(run):
        MODELDATA[inp] = recovershelve(INPUTS[inp]['ModelName'],[],recoverkeys)

#-------------------------------------------------------------------------------
# Plot Results
#-------------------------------------------------------------------------------
if True:

    #MomentPlot(MODELDATA,INPUTS,run1,'ParametricController_MomentSummary1',norm=True)
    #MomentPlot(MODELDATA,INPUTS,run2,'ParametricController_MomentSummary2',norm=True)
    MomentPlot(MODELDATA,INPUTS,run3,'ParametricController_MomentSummary3',norm=True)



    #now plot them.

    #make a 3d plotly showing moment diagrams from each modelset.

    #make a 3d plotly showing gider deflected shapes from each modelset

    #show reactions somehow. 

    #try to pull this plotly stuff from miser.