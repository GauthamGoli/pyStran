class Tsection():
    def __init__(self,BD,H,B,h,b,l,p):
        """ BD              :[startx,endx]
            Web Height      :H meters
            Flange width    :B
            Flange thickness:h
            Web Thickness   :b
            Lenght          :l
            density         :p kg/m3
            area            :A
            """
        
        self.start, self.end = float(BD[0]), float(BD[1])
        if (self.end-self.start)!= l:
            print "Wrong input, check again!"
            sys.exit()
        H,B,h,b,l,p=float(H),float(B),float(h),float(b),float(l),float(p)
        self.H, self.B, self.h, self.b, self.l, self.p =float(H), float(B), float(h), float(b), float(l), float(p)\
                                                         
        # cross sectional area
        self.A = B*h+H*b
        
        # y cog from bottom    
        self.ycog = ((H+h/2)*h*B+H**2*b/2)/self.A

        # Area moment of inertia Ixx
        self.I = b*H*(self.ycog-H/2)**2 + b*H**3/12 + h*B*(H + h/2 - self.ycog)**2 + h**3*B/12
        
        # Distance to Top surface from cog
        self.ytop = self.H+self.h-self.ycog
        
        # Distance to Bottom surface from cog
        self.ybottom = self.ycog
        
        if not isinstance(self,Tsection):
            self.ycog = self.H + self.h - self.ycog
            self.ytop = self.H + self.h - self.ycog
            self.ybottom = self.ycog
            
        # Mass
        self.M = self.A*self.l*self.p
        
        # Load related stuff all in KN/m
        selfload = self.M*9.81/(1000.0*self.l)
        liveload = 3000.0*9.81/(1000.0)
        self.load = selfload+liveload
        self.cost = self.A*self.l*2000.0/(0.3048**3)

class RTsection(Tsection):
    """ reverse of T section, takes input as T section and does calculations as reverse T section"""
    def __init__(self,BD,H,B,h,b,l,p):
        Tsection.__init__(self,BD,H,B,h,b,l,p)
        self.ycog = self.H + self.h - self.ycog
        self.ytop = self.H + self.h - self.ycog
        self.ybottom = self.ycog
