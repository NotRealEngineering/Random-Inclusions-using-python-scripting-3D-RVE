#Property of Not Real Engineering 
#Copyright 2020 Not Real Engineering - All Rights Reserved You may not use, 
#           distribute and modify this code without the written permission 
#           from Not Real Engineering.
############################################################################
##             Creating Random Inclusions                                 ##
############################################################################

from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random
from array import *
import math
import numpy
import os        # Operating system
import shutil    # copying or moving files

global rad;
rad = 2.0 # radius of inclusion

def partition(i,q):
    a=rad*sin(pi/6)
    b=rad*cos(pi/6)
## Creating Datum Planes 1
    dp1=mdb.models['Model-%d' %(q)].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=x_coordinate[i], principalPlane=YZPLANE)
    dp2=mdb.models['Model-%d' %(q)].parts['Part-1'].DatumPlaneByPrincipalPlane(offset=y_coordinate[i], principalPlane=XZPLANE)
	
## Creating Partition profile 2
    mdb.models['Model-%d' %(q)].ConstrainedSketch(gridSpacing=2.0, name='__profile__', sheetSize=25.0, transform=
        mdb.models['Model-%d' %(q)].parts['Part-1'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp1.id], sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, 5.0), ), sketchOrientation=TOP, origin=(x_coordinate[i],y_coordinate[i],z_coordinate[i])))
    mdb.models['Model-%d' %(q)].parts['Part-1'].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])
    c_1=mdb.models['Model-%d' %(q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(0,rad), point2=(rad,0))

    c_2=mdb.models['Model-%d' %(q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(rad,0), point2=(0,-rad))

    mdb.models['Model-%d' %(q)].sketches['__profile__'].Line(point1=(0,rad), point2=(0,-rad))
    mdb.models['Model-%d' %(q)].parts['Part-1'].PartitionCellBySketch(cells=
        mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((0.1, 0.1,0.1), )), sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'], 
        sketchOrientation=TOP, sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp1.id], sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, 5.0), ))
    
## Creating circle for giving path to sweep
    mdb.models['Model-%d' %(q)].ConstrainedSketch(gridSpacing=2.0, name='__profile__', sheetSize=25.0,
        transform=mdb.models['Model-%d' %(q)].parts['Part-1'].MakeSketchTransform(sketchPlane=mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp2.id],
        sketchPlaneSide=SIDE1,sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, 5.0), ), sketchOrientation=TOP, origin=(x_coordinate[i], y_coordinate[i], z_coordinate[i])))
    mdb.models['Model-%d' %(q)].parts['Part-1'].projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'])
    c_3=mdb.models['Model-%d' %(q)].sketches['__profile__'].CircleByCenterPerimeter(center=(0, 0), point1=(-rad, 0))

    mdb.models['Model-%d' %(q)].parts['Part-1'].PartitionCellBySketch(cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((0.2,0.2, 
        0.2), )), sketch=mdb.models['Model-%d' %(q)].sketches['__profile__'],sketchOrientation=TOP, sketchPlane=
        mdb.models['Model-%d' %(q)].parts['Part-1'].datums[dp2.id], sketchUpEdge=mdb.models['Model-%d' %(q)].parts['Part-1'].edges.findAt((0.0, 0.0, 5.0), ))
		
## Creating the spherical partition
    m= mdb.models['Model-%d' %(q)].parts['Part-1']
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((x_coordinate[i]+rad, y_coordinate[i], z_coordinate[i]),),cells=m.cells.findAt((0.2, 0.2, 0.2),),edges=(m.edges.findAt((x_coordinate[i], y_coordinate[i]-a,z_coordinate[i]+b), ),))
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((x_coordinate[i]+rad, y_coordinate[i], z_coordinate[i]),),cells=m.cells.findAt((0.2, 0.2, 0.2),),edges=(m.edges.findAt((x_coordinate[i], y_coordinate[i]+a,z_coordinate[i]+b), ),))	
    
dis=numpy.zeros(1000)

Max_iterations=4    # Set number of iterations
max_incl = 15      # set number of inclusions required

for q in range (1,Max_iterations):
    # LET'S CREATE MODEL
    mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-%d' %(q))    
    
    ## LETS CREATE MATRIX
    mdb.models['Model-%d' %(q)].ConstrainedSketch(name='__profile__', sheetSize=20.0)
    mdb.models['Model-%d' %(q)].sketches['__profile__'].sketchOptions.setValues(
        decimalPlaces=4)
    mdb.models['Model-%d' %(q)].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
        point2=(20.0, 20.0))
    mdb.models['Model-%d' %(q)].Part(dimensionality=THREE_D, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-%d' %(q)].parts['Part-1'].BaseSolidExtrude(depth=20.0, sketch=
        mdb.models['Model-%d' %(q)].sketches['__profile__'])
    del mdb.models['Model-%d' %(q)].sketches['__profile__']        

    num_incl = 0
    x_coordinate = []
    y_coordinate = []
    z_coordinate = []

    while (num_incl < max_incl):
        random_x=random.uniform(3.3, 16.7)
        random_y=random.uniform(3.3, 16.7)
        random_z=random.uniform(3.3, 16.7)

        isPointIntersecting = False
        for j in range (0,len(x_coordinate)):
    
    
            dis[j]=sqrt((random_x-x_coordinate[j])**2+(random_y-y_coordinate[j])**2+(random_z-z_coordinate[j])**2)

                
            if dis[j] < (2.2*rad):

                isPointIntersecting = True
                break

        if (isPointIntersecting == False):
            x_coordinate.append(random_x)
            y_coordinate.append(random_y)
            z_coordinate.append(random_z)
            num_incl = num_incl + 1        
    
    for i in range (num_incl):
        partition(i,q)

    # LET'S CREATE MATERIAL-1 (MATRIX POLYMER)
    mdb.models['Model-%d' %(q)].Material(name='Matrix')
    mdb.models['Model-%d' %(q)].materials['Matrix'].Elastic(table=
        ((1e2, 0.47), ))
    
    # LET'S CREATE MATERIAL-2 (ELASTIC INCLUSION)
    mdb.models['Model-%d' %(q)].Material(name='Elastic')
    mdb.models['Model-%d' %(q)].materials['Elastic'].Elastic(table=
        ((1e3, 0.35), ))
        
    # LET'S CREATE SECTIONS    
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='Matrix', name='Matrix', 
        thickness=None)
    mdb.models['Model-%d' %(q)].HomogeneousSolidSection(material='Elastic', name='Inclusion', 
        thickness=None)
        
    # LET'S ASSIGN SECTIONS
    mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((0.1, 0.1, 0.1), ), )), sectionName='Matrix',
        thicknessAssignment=FROM_SECTION)
    
    for i in range (num_incl):
        mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i], y_coordinate[i]-0.2*rad, z_coordinate[i]+0.2*rad), ), )),
            sectionName='Inclusion', thicknessAssignment=FROM_SECTION)
        mdb.models['Model-%d' %(q)].parts['Part-1'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            cells=mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i], y_coordinate[i]+0.2*rad, z_coordinate[i]+0.2*rad), ), )),
            sectionName='Inclusion', thicknessAssignment=FROM_SECTION)

    # LET'S CREATE INSTANCE
    mdb.models['Model-%d' %(q)].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-%d' %(q)].rootAssembly.Instance(dependent=ON, name='Part-1-1', 
        part=mdb.models['Model-%d' %(q)].parts['Part-1'])

    # LET'S CREATE STEP
    mdb.models['Model-%d' %(q)].StaticStep(initialInc=0.01, maxInc=0.1, maxNumInc=10000, 
        minInc=1e-12, name='Step-1', previous='Initial')
            
    # LET'S CREATE BOUNDARY CONDITIONS
    mdb.models['Model-%d' %(q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(
        faces=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.findAt(
        ((0.0, 1.0, 5.0), ), )), u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, 
        ur2=UNSET, ur3=UNSET)
    	
    mdb.models['Model-%d' %(q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(
        faces=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.findAt(
        ((1.0, 0.0, 5.0), ), )), u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, 
        ur2=UNSET, ur3=UNSET)				
    
    mdb.models['Model-%d' %(q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-3', region=Region(
        faces=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.findAt(
        ((5.0, 5.0, 0.0), ), )), u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, 
        ur2=UNSET, ur3=UNSET)

    mdb.models['Model-%d' %(q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-4', region=Region(
        faces=mdb.models['Model-%d' %(q)].rootAssembly.instances['Part-1-1'].faces.findAt(
        ((20.0, 5.0, 5.0), ), )), u1=0.1, u2=UNSET, u3=UNSET, ur1=UNSET, 
        ur2=UNSET, ur3=UNSET)

    # LET'S SEED THE PART 
    mdb.models['Model-%d' %(q)].parts['Part-1'].seedPart(deviationFactor=0.1, 
        minSizeFactor=0.1, size=1.0)
        
    # LET'S SET ELEMENT TYPE
    mdb.models['Model-%d' %(q)].parts['Part-1'].setElementType(elemTypes=(ElemType(
        elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
        elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
        regions=(mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((0.1,0.1,0.1), ), ), ))
        
    for i in range (num_incl):
        mdb.models['Model-%d' %(q)].parts['Part-1'].setElementType(elemTypes=(ElemType(elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
            elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
            regions=(mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i],y_coordinate[i]-0.2*rad,z_coordinate[i]+0.2*rad), ), ((x_coordinate[i],y_coordinate[i]+0.2*rad,z_coordinate[i]+0.2*rad), ), ), ))
            
    mdb.models['Model-%d' %(q)].parts['Part-1'].setMeshControls(elemShape=TET, regions=
        mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((0.1, 0.1, 0.1), ), ), technique=FREE)
    
    for i in range (num_incl):
        mdb.models['Model-%d' %(q)].parts['Part-1'].setMeshControls(elemShape=TET, regions=
    	    mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i],y_coordinate[i]-0.2*rad,z_coordinate[i]+0.2*rad), ), ((x_coordinate[i],y_coordinate[i]+0.2*rad,z_coordinate[i]+0.2*rad), ), ), technique=FREE)
        mdb.models['Model-%d' %(q)].parts['Part-1'].setMeshControls(elemShape=TET, regions=
    	    mdb.models['Model-%d' %(q)].parts['Part-1'].cells.findAt(((x_coordinate[i],y_coordinate[i]-1.05*rad,z_coordinate[i]), ), ((x_coordinate[i],y_coordinate[i]+1.05*rad,z_coordinate[i]), ), ), technique=FREE)                

    # LET'S GENERATE MESH
    mdb.models['Model-%d' %(q)].parts['Part-1'].generateMesh()

    #LET'S CREATE JOBS 
    mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
        explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
        memory=90, memoryUnits=PERCENTAGE, model='Model-%d' %(q), modelPrint=OFF, 
        multiprocessingMode=DEFAULT, name='Job-%d' %(q) , nodalOutputPrecision=SINGLE, 
        numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
        waitHours=0, waitMinutes=0)
    #mdb.jobs['Job-%d-%d' %(w,q)].writeInput()
    #mdb.jobs['Job-%d-%d' %(w,q) ].submit(consistencyChecking=OFF)    
    #mdb.jobs['Job-%d-%d' %(w,q) ].waitForCompletion()

    
#Property of Not Real Engineering 
