#from simplySupportedBeam import Beam,pointload,reaction,udl,uvl,roller,continuousBeam,hinge,moment,SimplySupportedBeam
from math import *
import numpy as np

class Node():
    def __init__(self,name,position,typ):
        self.id = name
        self.nodex = position[0]
        self.nodey = position[1]
        self.type = typ    #hinge,fixed joint,fixed end,free end,hinge joint
        self.counterNodes = []
        self.nodalforce = []  # can be a moment as well.
        self.constraints = {}
        self.displacements = {}
        self.d = {}  # holds selttlements and stuff?
        self.code_numbers = []  # keep track of DOF code numbers
        self.load_numbers = [] # Load DOF numbers
        
    def counternodes(self,counterNodes):
        #counterNodes = [2,5,3] list of counternodes
        self.counterNodes = counterNodes
        
    def displacement(self,number,direction,mag = 0):
        #specify displacement DOF and known displacements(if any)
        self.displacements[direction] = number
        self.code_numbers.append(number)
        if mag:
            self.d[number] = mag 
        
    def constraint(self,number,direction,mag = 0):
        self.constraints[direction] = number
        self.code_numbers.append(number)
        if mag:
            self.d[number] = mag

    def load(self,number,mag):
        loads[number] = mag
        self.load_numbers.append(number)

    

class Member():
    def __init__(self,node_N,node_F,A=1,E=1,I=1):
        self.near_end = node_N
        self.far_end = node_F
        self.A = A
        self.E = E
        self.I = I
        xn,yn,xf,yf=self.near_end.nodex,self.near_end.nodey,self.far_end.nodex,self.far_end.nodey
        self.length = sqrt((xf-xn)**2+(yf-yn)**2)
        self.lx = (xf-xn)/self.length
        self.ly = (yf-yn)/self.length
        self.force = None
        self.T = np.array([ [self.lx,self.ly,0,0,0,0],[-self.ly,self.lx,0,0,0,0],[0,0,1,0,0,0],[0,0,0,self.lx,self.ly,0],[0,0,0,-self.ly,self.lx,0],[0,0,0,0,0,1] ])
        self.Tt = self.T.T
        #setting up member stifness matrix
        k11 = self.A*self.E/self.length
        k22 = 12*self.E*self.I/(self.length)**3
        k33 = 4*self.E*self.I/self.length
        k32 = 6*self.E*self.I/(self.length)**2
        
        self.k1 = np.array([ [k11,0,0,-k11,0,0],[0,k22,k32,0,-k22,k32],[0,k32,k33,0,-k32,k33/2],[-k11,0,0,k11,0,0],[0,-k22,-k32,0,k22,-k32],[0,k32,k33/2,0,-k32,k33] ])  # member stiffness matrix
        #print len(self.Tt),len(self.k1),len(self.T)
        self.k = (self.Tt.dot(self.k1)).dot(self.T)
        displ1 = node_N.displacements.copy()
        displ2 = node_F.displacements.copy()
        displ1.update(node_N.constraints)
        displ2.update(node_F.constraints)
        self.nx,self.ny,self.nz,self.fx,self.fy,self.fz = displ1['x'],displ1['y'],displ1['z'],displ2['x'],displ2['y'],displ2['z']
        mapping = {0:self.nx,1:self.ny,2:self.nz,3:self.fx,4:self.fy,5:self.fz}
        mapping2 = mapping.copy()
        for point in mapping:
            for point2 in mapping2:
                K[mapping[point]-1,mapping2[point2]-1] += self.k[point,point2] 
    
    def load(self,load_type,parameters):
        # for member forces like an udl etc
        #LOAD_TYPE =1 means point load at any intermediate position              """ BUG CHECK LATER, ASSUMING NEAR END IS LEFT AND FAR END IS RIGHT FOR ALL CASES """"
        # parameters is a list [load magnitude, other parameters]
        if load_type=='pl':
            #parameters = [load,a,b]
            p,a,b = parameters[0],parameters[1],parameters[2]
            FEMn = p*b*b*a/(self.length)**2
            FEMf = p*b*a*a/(self.length)**2
            Vn = p*(1-a/self.length)-p*a*b*(a-b)/(self.length)**3
            Vf = p*a/self.length + p*a*b*(a-b)/(self.length)**3
        elif load_type=='udl':
            #parameters = [w]  UDL throughout
            w = parameters[0]
            FEMn = w*(self.length)**2/12.0
            FEMf = w*(self.length)**2/12.0
            Vn = w*self.length/2.0
            Vf = Vn
        elif load_type==3:
            #paramenters = [w] UDL for left half
            w = parameters[0]
            FEMn = 11*w*(self.length)**2/192.0
            FEMf = 5*w*(self.length)**2/192.0

            ##### SOMEONE FILL THIS SPACE #####
            
        elif load_type==4:
            # parameters = [w] left triangular load
            w = parameters[0]
            FEMn = w*(self.length)**2/20.0
            FEMf = w*(self.length)**2/30.0

            ### SOMEONE FILL THSI TOO ###

            
        elif load_type=='tr':
            # parameters = [w] central triangular load
            w = parameters[0]
            FEMn = 5*w*(self.length)**2/96.0
            FEMf = 5*w*(self.length)**2/96.0
            Vn = w*self.length/4.0
            Vf = Vn
            
        if 'z' in self.near_end.displacements:
            loads[self.near_end.displacements['z']] = -FEMn if self.near_end.displacements['z'] not in loads else loads[self.near_end.displacements['z']]-FEMn
            if self.near_end.displacements['z'] not in self.near_end.load_numbers:
                self.near_end.load_numbers.append(self.near_end.displacements['z'])
        if 'z' in self.far_end.displacements:
            loads[self.far_end.displacements['z']] = +FEMf if self.far_end.displacements['z'] not in loads else loads[self.far_end.displacements['z']]+FEMf
            if self.far_end.displacements['z'] not in self.far_end.load_numbers:
                self.far_end.load_numbers.append(self.far_end.displacements['z'])
        if 'y' in self.near_end.displacements:
            loads[self.near_end.displacements['y']] = -Vn*self.lx if self.near_end.displacements['y'] not in loads else loads[self.near_end.displacements['y']]-Vn*self.lx
            if self.near_end.displacements['y'] not in self.near_end.load_numbers:
                self.near_end.load_numbers.append(self.near_end.displacements['y'])
        if 'y' in self.far_end.displacements:
            loads[self.far_end.displacements['y']] = -Vn*self.lx if self.far_end.displacements['y'] not in loads else loads[self.far_end.displacements['y']]-Vn*self.lx
            if self.far_end.displacements['y'] not in self.far_end.load_numbers:
                self.far_end.load_numbers.append(self.far_end.displacements['y'])
        if 'x' in self.near_end.displacements:
            loads[self.near_end.displacements['x']] = Vn*self.ly if self.near_end.displacements['x'] not in loads else loads[self.near_end.displacements['x']]+Vn*self.ly
            if self.near_end.displacements['x'] not in self.near_end.load_numbers:
                self.near_end.load_numbers.append(self.near_end.displacements['x'])
        if 'x' in self.far_end.displacements:
            loads[self.far_end.displacements['x']] = Vn*self.ly if self.far_end.displacements['x'] not in loads else loads[self.far_end.displacements['x']]+Vn*self.ly
            if self.far_end.displacements['x'] not in self.far_end.load_numbers:
                self.far_end.load_numbers.append(self.far_end.displacements['x'])                                                
            
class Frame():
    def __init__(self,nodes,members):
        self.nodes = nodes
        self.members = members
        
    def solve(self):
        # call this after applying all the forces to set others to zero
        for node in self.nodes:
            if (node.type == 'joint' or node.type == 'fixedjoint') and len(node.load_numbers)!=len(node.code_numbers):
                for num in [n for n in node.code_numbers if n not in node.load_numbers]:
                    loads[num] = 0
            elif node.type == 'roller':
                loads[sorted(node.code_numbers)[0]] = 0
                loads[sorted(node.code_numbers)[1]] = 0   # PROBABLE BUG POINT DUE TO HARD CODING, CHECK LATER
        
        load = [[loads[l]] for l in sorted(loads)]
        Qk = np.array(load)
        K11 = K[:len(load),:len(load)]
        K12 = K[:len(load),len(load):]
        Dk = np.zeros((len(K)-len(K11),1))
        # fill Dk with settlements if any
        for node in self.nodes:
            for number in node.d:
                Dk[number-len(K11)-1] = node.d[number]

        Du = np.linalg.solve(K11,Qk-K12.dot(Dk))
        K21 = K[len(load):,:len(load)]
        K22 = K[len(load):,len(load):]
        Qu = K21.dot(Du)+K22.dot(Dk)
        self.Q = np.concatenate((Qk,Qu),axis=0)
        self.D = np.concatenate((Du,Dk),axis=0)

        """for member in self.members:
            member.force = member.A*member.E/member.length*np.array([ [-member.lx,-member.ly,member.lx,member.ly] ]).dot(np.array([ self.D[member.nx-1],self.D[member.ny-1],self.D[member.fx-1],self.D[member.fy-1] ]))     
            # print member.force                   """


#test code
nodes=[]
loads = {}
nodes.append(Node(1,[0,0],'fixed'))
nodes[0].constraint(4,'x')
nodes[0].constraint(6,'z')
nodes[0].constraint(5,'y')
nodes.append(Node(2,[240,180],'fixedjoint'))
nodes[1].displacement(2,'y')
nodes[1].displacement(1,'x')
nodes[1].displacement(3,'z')
nodes.append(Node(3,[480,180],'fixed'))
nodes[2].constraint(9,'z')
nodes[2].constraint(8,'y')
nodes[2].constraint(7,'x')
#nodes[3].load(3,2)
#nodes[3].load(4,-4)

dim=0
for node in nodes:
    dim += len(node.constraints)+len(node.displacements)
    
K = np.zeros((dim,dim))  #Global stiffness matrix (EMPTY) 

members=[]
members.append(Member(nodes[0],nodes[1],A=12,E=29000,I=600))
members.append(Member(nodes[1],nodes[2],A=12,E=29000,I=600))
members[1].load('udl',[0.25])
#members.append(Member(nodes[1],nodes[3]))

main = Frame(nodes,members)
main.solve()
