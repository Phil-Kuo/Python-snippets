from math import *
import numpy as np

class larCalc(object):
    def __init__(self, **kwargs):
        """ 
        Args
            self.tEnd: iteration time.
            self.Nt: # of time counts.
            self.tm: work time of motor.
        """
        # self.q = 0.0
        # self.r = 0.0
        # self.Sth = 0.0
        # self.Sps = 0.0
        # self.theta = 0.0
        # self.psi = 0.0
        self.tEnd = 200
        self.Nt = 1000
        self.dt = self.tEnd/self.Nt # simulation time interval
        self.tm = 30
        self.flight_time = np.arange(0, self.tEnd, self.Nt)

    def cmptThrust(self, P0):
        """ 
        compute Thrust.
        Args:
            P0: remain constant during motor working.
        return:
            P: shape(1, self.Nt)
        """
        P = []
        for t in self.flight_time:
            p = P0 if (t<self.tm) else 0
            P.append(p)        
        return np.ndarray(P)

    @staticmethod
    def cmptDrag(Vm, q, r):
        """
        compute drag force.
        Args:
            Vm: velocity of missile.
            q: pitch rate
            r: yaw rate
        Returns:            
            D: the missile's drag force.
        """
        K1 = 0.009412 # K1 = 0.5*rho*s*Cd0
        K2 = 93850/(9.8**2) # K2 = 2*K*w**2/(rho*s)
        D = K1*Vm**2 + K2 * (q**2 + r**2)/Vm**2
        return D
    
    @staticmethod
    def dynamicEqs(P, D, m, g0=9.8, Vm, theta, Sth, Sps, q, r, Cmt=0.25):
        """
        导弹的动力学方程组
        Args:
            P: Thrust.
            D: the missile's drag force.            
            m: mass of missile.
            g0: Gravity acceleration.
            Vm: velocity of missile
            theta: pitch angle.
            Sth: pitch signals.
            Sps: yaw signals
            q: pitch velocity.
            r: yaw velocity.
            Cmt: missile time constant, =0.25s.
            Ne: # of equations.
        Returns:
            y: Equations for integration.
            index           Value
            0     Vmdot: acceleration of missile 
            1     thetadot: first derivative of theta vesus time 
            2     psidot: first derivative of psi versus time.
            3     qdot: pitch rate
            4     rdot: yaw rate
        """
        y = np.zeros((5, 1)) # 初始化

        y[0] = (P - D)/m -g0*sin(theta)
        y[1] = (q - cos(theta))/Vm
        y[2] = r/(Vm*cos(theta))
        y[3] = (Sth - q)/Cmt
        y[4] = (Sps - r)/Cmt
        return y            

    def odeRK4(self, y0, P, D, m, g0=9.8, Vm, theta, Sth, Sps, q, r, Cmt=0.25):
        """
        利用四阶龙格库塔法积分运算
        Args:           
           y0: 初始值
        returns:
            y: y[0,:], y[1,:], y[2,:], y[3,:], y[4,:]分别代表Vm、theta、psi、q、r
        """
        y = np.zeros((5, self.Nt))
        y[:,0] = y0 # 初值        
        
        for i in range(self.Nt)
            k1 = self.dynamicEqs(y[:,i])
            k2 = self.dynamicEqs(y[:,i] + self.dt * k1/2)
            k3 = self.dynamicEqs(y[:,i] + self.dt * k2/2)
            k4 = self.dynamicEqs(y[:,i] + self.dt * k3)
            y[:,i+1] = y[:,i] + self.dt/6 * (k1 + 2*k2 +2*k3 +k4)            
        return y

    def odeMotions(self, V, theta, psi):
        """
        compute position vector by integration.
        Args: 
            V: shape(1, self.Nt), velocity.
            theta: shape(1, self.Nt), pitch angle.
            psi: shape(1, self.Nt), yaw angle.
        returns:
            R: shape(3, self.Nt), positon vector, [X; Y; Z]
        """        
        
        R = np.zeros((3, self.Nt))
        
        dx = V*cos(theta)*cos(psi)*self.dt
        dy = V*cos(theta)*sin(psi)*self.dt
        dz = V*sin(theta)*self.dt

        R[0, :] += dx
        R[1, :] += dy
        R[2, :] += dz

        R = np.cumsum(R, axis=0)
        return R

    def cmptMass(self, mMax, mMin):
        """ 
        compute mass, considering mass lost due to oxidation.  
        Args:
            mMax: maximun mass.
            mMin: mass after burn-out.           
        return:
            m: mass, shape(1, self.Nt)
        """
        k = (mMax- mMin)/self.tm
        m = []
        for (i, t) in zip(np.range(self.Nt), self.flight_time):
            m[i] = (mMax- k*t) if (t<self.tm) else mMin
        return np.ndarray(m)

    @staticmethod
    def cmptDist(Rm, Rt):
        """
        Args:
            Rm: shape(3, self.Nt), [xm; ym; zm], position vector of missile.
            Rt: shape(3, self.Nt), [xt; yt; zt], position vector of target.
        return:
            deltaDist: shape(1, self.Nt), distance of missile and target.
            delta: shape(3, self.Nt), [deltaX; deltaY; deltaZ], difference of target and missile.
        """
        delta = Rt - Rm

        deltaDist = np.sqrt(np.sum((Rm - Rt)**2, axis=0))
        return deltaDist, delta

    def cmptClsVel(self, deltaDist, delta):
        """ 
        compute closing velocity
        Args:
            deltaDist: shape(1, self.Nt), distance of missile and target.
            delta: shape(3, self.Nt), [deltaX; deltaY; deltaZ], difference of target and missile.
        return:
            Vc: shape(1, self.Nt-1), closing velocity.
            CR: shape(3, self.Nt-1), closing rate, [CRx; CRy; CRz], compute closing rate in X, Y, Z positions.
        """       
        CR = -np.diff(delta, axis=1)/self.dt

        Vc = np.diff(deltaDist)/self.dt
        return Vc, CR

    def cmptDesireVel(self, CR, delta, deltaDist, theta, psi):
        """
        compute desire velocities in position X, Y, Z.
        Args:            
            CR: shape(3, self.Nt-1), closing rate, [CRx; CRy; CRz], compute closing rate in X, Y, Z positions.           
            delta: shape(3, self.Nt), [deltaX; deltaY; deltaZ], difference of target and missile.
            deltaDist: shape(1, self.Nt), distance of missile and target.
            theta: shape(1, self.Nt), pitch angle.
            psi: shape(1, self.Nt), yaw angle.
        return:
            Vdesire: shape(3, self.Nt-1), desired velocities in the body axis.
            qDesire: shape(1, self.Nt-1), desired pitch rate.
            rDesire: shape(1, self.Nt-1), desired yaw rate.
        """        
        Vdesire = []
        for i in range(self.Nt-1):
            Vdesire[i] = np.cross(delta[:,i], CR[:,i])/deltaDist[i]**2 
        Vdesire = np.ndarray(Vdesire)
        qDesire = -np.multiply(sin(psi), Vdesire[0,:]) + np.multiply(cos(psi), Vdesire[1,:])
        rDesire = np.multiply(cos(theta), Vdesire[1,:]) + np.multiply(sin(theta), np.multiply(cos(psi, Vdesire[0,:])+np.multiply(sin(psi,Vdesire[1,:]))))
        return Vdesire, qDesire, rDesire

    @staticmethod
    def cmptsignals(Cnav=4, Vc, qDesire, rDesire):
        """ 
        compute the pitch and yaw signals.
        Args:
            Cnav: guidance constants, equals 4.
            Vc: shape(1, self.Nt-1), closing velocity.
            qDesire: shape(1, self.Nt-1), desired pitch rate.
            rDesire: shape(1, self.Nt-1), desired yaw rate.
        return:
            Sth: shape(1, self.Nt-1), the pitch(theta) signal.
            Sps: shape(1, self.Nt-1), the yaw(psi) signal.
        """
        Sth = Cnav * np.multiply(Vc, qDesire)
        Sps = Cnav * np.multiply(Vc, rDesire)
        return Sth, Sps


    def simulation(self, parameter_list):
        
        # positions of missile and target.
        self.odeMotions(Vm, theta, psi)
        self.odeMotions(Vt, thetat, psit)

        # calculate pitch and yaw signals for integration.
        deltaDist, delta = self.cmptDist(Rm, Rt)
        Vc, CR = self.cmptClsVel(deltaDist, delta)
        Vdesire, qDesire, rDesire = self.cmptDesireVel(CR, delta, deltaDist, theta, psi)
        Sth, Sps = self.cmptsignals(Cnav=4, Vc, qDesire, rDesire)

        P = self.cmptThrust(P0)
        D = self.cmptDrag(Vm, q, r)
        m = self.cmptMass(mMax, mMin, tm)

        y0 = [200; 0; 0; 0; 0]
        y = self.odeRK4(y0)