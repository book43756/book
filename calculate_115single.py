import numpy as np
from icecream import ic

class single115:
    def __init__ (self,Xp,Yp,Size):
        self.Xp = Xp
        self.Yp = Yp
        self.Size = Size
        # * กำหนดตำแหน่งเฟส A,B,C
        self.Xa = 0; 
        self.Ya = 15.1
        self.Xb = 4  
        self.Yb = 15.1
        self.Xc = 4  
        self.Yc = 12.6
        # * กำหนดแรงดันและมุมเฟส A,B,C
        self.r_a = 115
        self.theta_a = 0
        self.r_b = 115
        self.theta_b = -120
        self.r_c = 115
        self.theta_c = 120
        self.E0=8.854*pow(10,-12)
        
    def pol2cart(self,r,theta):
        z = r * np.exp(1j * theta)
        x, y = z.real, z.imag
        return x, y
    def cart2pol(self,x, y):
        z = x + y * 1j
        r,theta = np.abs(z), np.angle(z)*57.2958
        #! 1rad = 57.2958 deg 
        return r,theta
    def cal(self):
        self.Xva,self.Yva = self.pol2cart(self.r_a,np.radians(self.theta_a))
        Xvb,Yvb = self.pol2cart(self.r_b,np.radians(self.theta_b))
        Xvc,Yvc = self.pol2cart(self.r_c,np.radians(self.theta_c))
        V=np.array([[complex(Xva,Yva)/np.sqrt(3)],[complex(Xvb,Yvb)/np.sqrt(3)],[complex(Xvc,Yvc)/np.sqrt(3)]])

Test1 = single115(1,2,3)
Test1.cal()
