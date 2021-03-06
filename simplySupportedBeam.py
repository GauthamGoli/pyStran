# Able to solve indeterminate continuous BEAMS. Find all the reactions

import matplotlib.pyplot as plt
from operator import itemgetter
import numpy,sys
from sections import Tsection

class Beam():
    """
    base class for a beam
    """
    def __init__(self,b_length):
        self.length = b_length
        
class pointload():
    def __init__(self,x,f_load):
        self.startx = float(x)
        self.load = float(f_load)
    def getLoad(self,x=None):
        return self.load
    def moment(self,x):
        x = float(x)
        if x<self.startx:
            return self.getLoad()*(self.startx-x)
        if x>self.startx:
            return self.getLoad()*(x-self.startx)
        else:
            return 0.0
    

class reaction(pointload):
    """ reactions are same as pointloads, but they are presumed to be upward and hence their computed value is positive
    but while calculating bmd the moment should be negative, hence the getLoad() is overridden"""
    def getLoad(self,x=None):
        return -self.load
    

class udl():
    """ udl class """
    def __init__(self,x_start,x_end,f_load):
        self.startx = float(x_start)
        self.endx = float(x_end)
        self.span = self.endx - self.startx
        self.load = float(f_load)
    def moment(self,x,dirn=None):
        """ Returns moment of udl about a point """
        x = float(x)
        if x>= self.endx:
            return self.getLoad()*(x-self.endx+self.span/2)
        elif x>self.startx and x<self.endx:
            return self.getLoad(x)*(x-self.startx)/2
        else:
            return 0
    def getLoad(self,x=None):
        """ returns magnitude of load"""
        if x and x<=self.endx:
            return self.load*(x-self.startx)
        return (self.span)*self.load
        
class uvl():
  #  uniformly varying load
    def __init__(self,config):
      #  """ config : [[x,load_magnitude],[x2,load_magnitude2],....]"""
        self.config = sorted(config,key=itemgetter(0))
        for jo in self.config:
            jo[0],jo[1] = float(jo[0]),float(jo[1])
        self.startx = float(self.config[0][0])
        self.endx = float(self.config[-1][0])

    def moment(self,x):
     # Returns moment towards left of the point x
        if x == self.startx:
            return 0
       # print "self.getload("+str(x)+") = "+ str(self.getLoad(x))+ '*' + "self.getcenroid() " + str(self.getCentroid(x))
        return self.getLoad(x)*self.getCentroid(x)


    def getInterHeight(self,seg_0,seg_1,width,widthDiff,heightDiff):
        if seg_1[1]>seg_0[1]:
            interheight = seg_0[1] + width * heightDiff / widthDiff
        else:
            interheight = seg_1[1] + ( widthDiff - width )*heightDiff / widthDiff
        return interheight

    
    def getLoad(self,x=None):
        # returns load upto left of x , default argument when we just want the total load of uvl
        x = self.endx if x == None else x
        segment_0,area = None, 0.0
        if x == self.startx:
            return 0.0    
        for segment in self.config:
            if segment[0]<=x and not segment_0:
                segment_0 = segment
                continue
            if segment[0]>x and segment_0[0]<x:
                width = float(x-segment_0[0])
                widthDiff = float(segment[0] - segment_0[0])
                heightDiff = abs(segment[1]-segment_0[1])
                interheight = self.getInterHeight(segment_0,segment,width,widthDiff,heightDiff)
                averageheight = ( segment_0[1] + interheight ) / 2
                area += averageheight*width
                break
            if segment!=segment_0:
                width = float(segment[0] - segment_0[0])
                averageheight = (segment[1]+segment_0[1])/2.0
                area += averageheight * width
                segment_0 = segment
        #print "area = " + str(area)
        return area

    
        

    def getCentroid(self,x):
        # distance of centroid of uvl configuration left of x, from x
        centroid,avg,seg_0 = 0.0,0,None
        for segment in self.config:
            if not seg_0:
                seg_0 = segment
                continue
            elif segment[0]<=x:
                h1, h2, w = float(seg_0[1]), float(segment[1]), float((segment[0]-seg_0[0]))
                centroid = x - ((h1 + 2*h2)*(w)/(3.0*(h1+h2)) + seg_0[0])
                break
            elif segment[0] > x:
                h1,h2, w = seg_0[1],self.getInterHeight(seg_0,segment,(x-seg_0[0]),(segment[0]-seg_0[0]),abs((segment[1]-seg_0[1]))),(x-seg_0[0])
                centroid = x - ((h1 + 2*h2)*(w)/(3*(h1+h2)) + seg_0[0])
                break
            
        #centroid = x - float(centroid/avg)
        #print "centroid= " + str(centroid)
        return centroid
                        
                    
class moment():
    """ moment class : anticlockwise is positive"""
    def __init__(self,x,torque):
        self.startx = float(x)
        self.torque = float(torque)
    def moment(self,x,dirn=None):
        if (x>=self.startx and self.torque<0):
            return self.torque
        elif x>=self.startx and self.torque>0:
            return -self.torque
    
            
class SimplySupportedBeam():
    """
    class for Simplysupported beam
    """
    
    def __init__(self,b_length,b_e=None,b_i=None,i_Type = None):
        #i_Type can be either 'linear' or 'constant' or 'poly'
        self.length = float(b_length)
        #print 'type' + str(type(b_e))
        self.pointloads,self.udl,self.moments,self.reactions,self.uvl,self.e = [],[],[],[],[],b_e
        self.iType = i_Type
        self.i = b_i
        self.My,self.Sy,self.x_samples = [],[],[]
        self.connections,self.supports = [],[]
        self.parts = 10000
        #print "Enter point loads as [ [location,load in kN],.... ] and call the method applyPointLoads() "
        #print "Enter udl as [[start,end,load in kN/m ]....] and call the method applyUDL()"
        #print "Enter uvl as [[x,load_magnitude],[x2,load_magnitude2],....] and call the method applyUVL()"
        # UNITS OF E: GPa, UNITS OF I: mm4
    def applyPointLoads(self, b_loads):
        """
        b_loads = [ [location,load in kN],.... ]
        """
        for entry in b_loads:
            self.pointloads.append(pointload(entry[0],entry[1]))
    def applyUDL(self, b_udl):
        """
        b_udl = [[ start , end , load in kN ]....]
        """
        for entry in b_udl:
            self.udl.append(udl(entry[0],entry[1],entry[2]))

    def applyUVL(self, b_uvl):
        """
        b_uvl = [[x,load_magnitude],[x2,load_magnitude2],....]
        """
        self.uvl.append(uvl(b_uvl))
        
    def applyMoment(self, b_moments):
        """
        b_moments = [[x,moment(anticlock +)],....]
        """
        for entry in b_moments:
            x, torque = entry[0], entry[1]
            self.moments.append(moment(x,torque))
    def getEI(self,x):
        """ returns EI Value at a particular x. UNITS I: mm**4, E: GPa """
        if self.iType == 'linear':
            """DS of I( if varying linearly between x1 and x2) : [[x1,I1],[x2,I2]...] """
            I1, I2, x1, x2 = float(self.i[0][1]),float(self.i[1][1]),float(self.i[0][0]),float(self.i[1][0])
            slope = (I2-I1)/(x2-x1)
            Ix = (I1 + (slope)*(x-x1))*10**-6
            return Ix*(self.e)
        elif self.iType == 'constant':
            """ if it is const DS OF I : [[x1,x2,I1],[x3,x4,I2]...]"""
            for limit in self.i:
                iStart,iEnd,Ix = limit[0],limit[1],limit[2]*10**-6
                if iStart<=x and iEnd>x or self.length == x:
                    return Ix*(self.e)
                continue
        elif self.iType =='poly':
            """ if it is parabolic: [x=[array of x values],y=[array of y values(depth? may be or I values)]] example: b_i=[[0,5,10],[-3,0,3]]  """
            p = numpy.poly1d(numpy.polyfit(x,[i*10**-12 for i in y],2))
            Ix = p(x)
            return Ix*(self.e)*10**-6
            
            
                               
    def findreactions(self):
        """ prints reactions of the ssBeam"""
        rA,rB = 0,0
        for load in self.pointloads+self.udl+self.moments+self.uvl:
            rA += load.moment(self.length)
        #print rA
        rA = rA/self.length
        for load in self.pointloads+self.udl+self.uvl:
            rB += load.getLoad()
        rB -= rA
        self.reactions.append(reaction(0,rA))
        self.reactions.append(reaction(self.length,rB))
        print 'Reactions: rA(left end) = %s, rB(right end) = %s' % (rA,rB)
        
    def findBM(self,x):
        """ finds the Bending Moment at a point on beam """
        required,bm = [], 0 
        for force in self.pointloads+self.moments+self.reactions+self.udl+self.uvl:
            if force.startx <= x:
                bm += force.moment(x)
        #print "bm at %s = %s " %(x,bm)
        return bm
    def findSF(self,x):
        """ finds the shear force at a point on beam """
        sf = 0
        for force in self.pointloads+self.reactions + self.udl + self.uvl:
            if x==0:
                return 0
            if force.startx <= x:
                sf += force.getLoad(x)
        #print "sf at %s = %s " %(x,sf)
        return sf

    def calculations(self):
        x,y,x1,y1 = [], [], [], []
        #for i in range(0,int(self.length*1000)+5,5):
        # self.parts is number of parts, value also will be used in trapz integration - deflection().
        for i in numpy.linspace(0,self.length,num=self.parts,endpoint=True):
            x.append(i)
            y.append(-self.findSF(x[-1]))
            x1.append(i)
            y1.append(-self.findBM(x1[-1]))
        self.My,self.Sy,self.x_samples = y1,y,x
        
    def plotBMSFD(self):  
        f, axarr = plt.subplots(2, sharex=True)
        plotGraph(self.x_samples,self.Sy,'length','magnitude','ShearforceDiagram',axarr[0])
        plotGraph(self.x_samples,self.My,'length','magnitude','Bending Moment Diagram',axarr[1])
        plt.show()

    def plotBDD(self):
        """Plots flexural compression and tension along extreme fibers BDD:Boundary Diagrams :P"""
        f, axarr = plt.subplots(5, sharex=True)
        tension,compression,shear = [],[],[]
        for x,My in zip(self.x_samples,self.My):
            section = [sect for sect in sections if (sect.start==0 and x==0) or (sect.start<x and sect.end>=x)]
            if len(section)!=1:
                print x
                print "check the code dude!"
                sys.exit()
            section = section.pop()
            if My>=0: #sagging
                tension.append(My*section.ybottom*10**-3/section.I)
                compression.append(My*section.ytop*10**-3/section.I)
            elif My<0:
                compression.append(-My*section.ybottom*10**-3/section.I)
                tension.append(-My*section.ytop*10**-3/section.I)
        for x,Sy in zip(self.x_samples,self.Sy):
            section = [sect for sect in sections if (sect.start==0 and x==0) or (sect.start<x and sect.end>=x)]
            if len(section)>1:
                print "check code"
                sys.exit()
            section = section.pop()
            shear.append(abs(Sy)*10**-3/section.A)
        plotGraph(self.x_samples, tension,'','magnitude MPa','FlexuralTension',axarr[0])
        plotGraph(self.x_samples, compression,'','magnitude MPa','FlexuralCompression',axarr[1])
        plotGraph(self.x_samples, shear,'','magnitude MPa','ShearStress',axarr[2])
        plotGraph(self.x_samples,self.Sy,'','magnitude kN','ShearforceDiagram',axarr[3])
        plotGraph(self.x_samples,self.My,'length','magnitude kN','Bending Moment Diagram',axarr[4])
        plt.show()
        
    def deflection(self,x):
            """ finds  the deflection at a point using unit load method, My is an array of M values along beam"""
            if isinstance(self,continuousBeam):
                beam = continuousBeam(self.length)
                beam.connections,beam.supports = self.connections,self.supports
            else:
                beam = SimplySupportedBeam(self.length,self.e,self.i)
            beam.applyPointLoads([[x,1.0]])
            beam.findreactions()
            beam.calculations()
            if not self.reactions:
                self.findreactions()
                self.calculations()
            #beam.plotBMSFD()
            I = []
            for i in range(len(self.x_samples)):
                try:
                    I.append(beam.My[int(i)]*self.My[int(i)]/self.getEI(self.x_samples[i]))
                    #I.append(beam.My[i]*self.My[i])
                    #print 'my= %s My = %s ' %(beam.My[i],self.My[i])
                except Exception,e:
                    print 'EI Values not specified for x= %s, exiting...with following error : %s' %(self.x_samples[int(i)-1],e)
                    return
            defl = numpy.trapz(I,dx=self.length/self.parts)
            print "Deflection at x = %s = %s m" %(x,defl)
            return defl

        
def plotGraph(xdata,ydata,x_label,y_label,header,plto):
    plto.plot(xdata,ydata)
    plto.set_title(header)
    plto.set_xlabel(x_label)
    plto.set_ylabel(y_label)
    plto.plot([0,xdata[-1]],[0,0],color = "black",linewidth=4.0)
    plto.fill_between(xdata,0,ydata,color="none",hatch="///",edgecolor="b")


class pin():
    def __init__(self,x):
        self.x = x
        self.type = 'pin'
class hinge():
    def __init__(self,x):
        self.x = x
        self.type = 'hinge'
class roller():
    def __init__(self,x):
        self.x = x
        self.type = 'roller'

class continuousBeam(SimplySupportedBeam):
    def specifySupports(self,pinArray,rollerArray,hingeArray):
        """pinArray is an array of x coordinates for pin supports
           similarly for rollerArray and hingeArray """
        for xpos in sorted(pinArray):
            self.supports.append(pin(xpos))
        for xpos in sorted(hingeArray):
            self.connections.append(hinge(xpos))
        for xpos in sorted(rollerArray):
            self.supports.append(roller(xpos))

    def checkstability(self):
        """ checks stablility of the continuous beam (DONT USE THIS . Lost track of this and is completely wrong. rewrite this """
        hinges = [connection for connection in self.connections if connection.type=='hinge']
        for hinge in hinges:
            print hinge.x
        for hinge in hinges:
            flag = True # assuming the structure is unstable
            #right,left = rightsupports(hinge.x),leftsupports(hinge.x)
            right,left = [support for support in self.supports if support.x>hinge.x and support.type!='hinge'],[support for support in self.supports if support.x<hinge.x and support.type!='hinge']
            for support in right:
                print support.type
                print support.x
            for support in right:
                if support.type!='roller':
                    flag = False
            if flag==True: # the beam is unstable
                print 'The beam is unstable at x = %s' %(right[0].x)
            for support in left:
                if support.type!='roller':
                    flag = False
            if flag==True:
                print 'The beam is unstale at x = %s' %(left[0].x)
        
    def rightsupports(xpos):
        right = [support for support in self.supports if support.x>xpos and support.type!='hinge']
        #locations = sorted([support.x for support in right])
        return right
    
    def leftsupports(xpos):
        left = [support for support in self.supports if support.x<xpos and support.type!='hinge']
        return left
            
        
    def findreactions(self):
        totalload,a,b = 0,[],[]
        for load in self.pointloads+self.udl+self.uvl+self.reactions:
            totalload += load.getLoad()
        b.append([totalload])
        """ b matrix """
        for connection in self.connections:
            bm = 0
            if connection.type == 'hinge':
                for force in self.pointloads+self.moments+self.udl+self.uvl+self.reactions:
                    if force.startx <=connection.x:
                        bm += force.moment(connection.x)
                b.append([bm])
        b.append([self.findBM(self.length)])
        """ a matrix """
        """ if some reactions are already found (Indeterminate beam) exclude them from the equations"""
        excludedSupports = [reac.startx for reac in self.reactions]
        a.append([1 for support in self.supports if support.type!='hinge' and support.x not in excludedSupports])
        for connection in self.connections:
            nextarr = [connection.x-support.x for support in self.supports if connection.type == 'hinge' and connection.x>support.x and support.x not in excludedSupports] # next row matrix in a
            nextarr.extend([0]*(len(b)-len(nextarr))) # filling the remaining elements with zeros
            a.append(nextarr)
        nextarr = [self.length-support.x for support in self.supports if support.x not in excludedSupports]
        a.append(nextarr)
        if len(a) != len(b):
            print 'run run!'
        elif len(a)>len(a[0]):
            print 'The structure is unstable(by preliminary examination)! Re design the structure'
            sys.exit()
        elif len(a)<len(a[0]):
            print 'The structure is indeterminate! Try Solving it?(y/n)'
            res = raw_input('>')
            if res=='y':
                print '%d redundant supports may be there' %(len(a[0])-len(a))
                self.indeterminate()
                return
            else:
                sys.exit()
        X = numpy.linalg.solve(numpy.array(a),numpy.array(b))
        self.reactions.extend([reaction(support.x,reac[0]) for reac,support in zip(X,self.supports)])
        #print ['Reactions at x = %s = %s kN' %(support.x,reac[0]) for reac,support in zip(X,self.supports)]
        print ['Reactions at x = %s = %s kN' %(reac.startx,-reac.getLoad()) for reac in self.reactions]

    def indeterminate(self):
        print 'Enter the positions of redundant supports(x1,x2,x3...):'
        res = raw_input('>')
        red = sorted([float(e) for e in res.split(',')])
        self.supports = [support for support in self.supports if support.x not in red] #removing the redundant supports
        self.findreactions()
        self.calculations()
        b,a = [self.deflection(x) for x in red],[[] for e in red]  # b matrix , empty  a matrix
        #cloning the beam, except for loadings and redundant supports
        rBeam = continuousBeam(self.length,self.e,self.i,self.iType)
        rBeam.supports,rBeam.connections=self.supports,self.connections 
        for x1 in red:
            for x2 in red:
                rBeam.pointloads=[pointload(x2,1)]
                rBeam.findreactions()
                rBeam.calculations()
                a[red.index(x1)].append(rBeam.deflection(x1))
        X = numpy.linalg.solve(numpy.array(a),numpy.array(b))
        print ['Reactions at x = %s =%s kN' %(support,reac) for reac,support in zip(X,red)]
        self.reactions=[]
        for x,reac in zip(red,X):
            self.reactions.append(reaction(x,reac))
        self.findreactions()
            
    
#test code    case 0
"""beam1 = continuousBeam(b_length = 8,b_e = 9.0,b_i = [[0,4.0,14.6*10**9],[4.0,8.0,(14.6*2)*10**9]], i_Type = 'constant')
beam1.specifySupports(pinArray=[0.0],rollerArray=[8.0],hingeArray=[])
#print beam1.getEI(7)
beam1.applyUDL([[0,8,8]])
beam1.applyPointLoads([[4,20]])
#beam1.applyPointLoads([[7.5,8]])
beam1.findreactions()
beam1.calculations()
beam1.deflection(4.0)"""
#case1
"""
beam1 = continuousBeam(20)
beam1.specifySupports(pinArray=[0,5,15],rollerArray=[10.0,20],hingeArray=[12,17])
beam1.checkstability()
beam1.applyUDL([[0,4,10]])
beam1.applyPointLoads([[10,20]])
beam1.findreactions()
beam1.calculations()
beam1.plotBMSFD()"""

#case2
cost = 2000.0/(0.3048**3)
sections=[]
sections.append(Tsection([0,37],H=0.25,B=4,h=0.3,b=2,l=37,p=600))
sections.append(Tsection([37,120],H=0.4,B=4,h=0.4,b=2,l=83,p=600))

beam1 = continuousBeam(120,b_e=10,b_i=[[0,37.0,49500000256],[37,120,117333336064]],i_Type = 'constant')
beam1.specifySupports(pinArray=[0,37],rollerArray=[97],hingeArray=[37.0])
#beam1.checkstability()
beam1.applyUDL([[0,37,sections[0].load],[37,120,sections[1].load]])
#beam1.applyPointLoads([[17,0.5]])
beam1.findreactions()
beam1.calculations()
#beam1.plotBMSFD()

beam1.plotBDD()

# testcode        

#beam1 = SimplySupportedBeam(8,200,[[0,400*10**6],[8,400*10**6]])
#beam1.applyUDL([[4,8,10]])
#beam1.applyUDL([[12,16,2]])
#beam1.applyPointLoads([[2.0,5]])
#beam1.applyPointLoads([[6,6]])
#beam1.applyUVL([[0,5],[6,3]])
#beam1.applyUVL([[7,6],[14,0]])
#beam1.applyMoment([[0,5]])
#beam1.findreactions()
#beam1.calculations()
#beam1.plotBMSFD()
#x,y = [], []
"""for i in range(0,int(beam1.length*100)+5,5):
    x.append(i/100.0)
    y.append(-beam1.findSF(x[-1]))
x1,y1 = [],[]
for i in range(0,int(beam1.length*100)+5,5):
    x1.append(i/100.0)
    y1.append(-beam1.findBM(x1[-1]))
    #print "x="+str(x[-1])+": "+str(y[-1]) """
#beam1.deflection(4.0)
#print 'def test= ' + str(beam1.deflection(4.0))
#f, axarr = plt.subplots(2, sharex=True)  
#plotGraph(x,y,'length','magnitude','ShearforceDiagram',axarr[0])

#plotGraph(x1,y1,'length','magnitude','Bending Moment Diagram',axarr[1])

#plt.show()
#plt.savefig('graph1.png')
#from simplySupportedBeam import *
#import numpy

