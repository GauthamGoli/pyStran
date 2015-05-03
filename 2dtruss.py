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
    def __init__(self,node_N,node_F,A=1,E=1):
        self.near_end = node_N
        self.far_end = node_F
        self.A = A
        self.E = E
        xn,yn,xf,yf=self.near_end.nodex,self.near_end.nodey,self.far_end.nodex,self.far_end.nodey
        self.length = sqrt((xf-xn)**2+(yf-yn)**2)
        self.lx = (xf-xn)/self.length
        self.ly = (yf-yn)/self.length
        self.force = None
        self.T = np.array([ [self.lx,self.ly,0,0],[0,0,self.lx,self.ly] ])
        self.Tt = self.T.T
        self.k1 = self.A*self.E/self.length * np.array([ [1,-1],[-1,1] ])  # member stiffness matrix
        self.k = (self.Tt.dot(self.k1)).dot(self.T)
        displ1 = node_N.displacements.copy()
        displ2 = node_F.displacements.copy()
        displ1.update(node_N.constraints)
        displ2.update(node_F.constraints)
        self.nx,self.ny,self.fx,self.fy = displ1['x'],displ1['y'],displ2['x'],displ2['y']
        mapping = {0:self.nx,1:self.ny,2:self.fx,3:self.fy}
        mapping2 = mapping.copy()
        for point in mapping:
            for point2 in mapping2:
                K[mapping[point]-1,mapping2[point2]-1] += self.k[point,point2] 
    

    

class Truss():
    def __init__(self,nodes,members):
        self.nodes = nodes
        self.members = members
        
    def solve(self):
        # call this after applying all the forces to set others to zero
        for node in self.nodes:
            if node.type == 'joint' and len(node.load_numbers)!=len(node.code_numbers):
                for num in [n for n in node.code_numbers if n not in node.load_numbers]:
                    loads[num] = 0
            elif node.type == 'roller':
                loads[node.code_numbers[0]] = 0
            
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

        for member in self.members:
            member.force = member.A*member.E/member.length*np.array([ [-member.lx,-member.ly,member.lx,member.ly] ]).dot(np.array([ self.D[member.nx-1],self.D[member.ny-1],self.D[member.fx-1],self.D[member.fy-1] ]))     
            # print member.force                   


#test code
nodes=[]
loads = {}
nodes.append(Node(1,[4,0],'hinge'))
nodes[0].displacement(3,'x')
nodes[0].displacement(4,'y',-0.025)
nodes.append(Node(2,[4,3],'joint'))
nodes[1].displacement(2,'y')
nodes[1].constraint(1,'x')
nodes.append(Node(3,[0,0],'hinge'))
nodes[2].constraint(5,'x')
nodes[2].constraint(6,'y')
nodes.append(Node(4,[0,3],'hinge'))
nodes[3].displacement(7,'x')
nodes[3].displacement(8,'y')
#nodes[3].load(3,2)
#nodes[3].load(4,-4)

dim=0
for node in nodes:
    dim += len(node.constraints)+len(node.displacements)
    
K = np.zeros((dim,dim))  #Global stiffness matrix (EMPTY) 

members=[]
members.append(Member(nodes[0],nodes[1]))
members.append(Member(nodes[1],nodes[2]))
members.append(Member(nodes[1],nodes[3]))

main = Truss(nodes,members)
main.solve()
