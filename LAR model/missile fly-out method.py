from math import *
import numpy as np

class larCalc(object):
    def __init__(self, parameter_list):
        pass

    @staticmethod
    def cmptDrag(Vm, q, r):
        """
        计算阻力

        Args:
            Vm: velocity of missile.
            q: pitch rate
            r: yaw rate
        Returns:            
            D: the missile's drag force.
        """
        K1 = 0.5*rho*s*Cd0
        K2 = 2*K*w**2/(rho*s)
        D = K1*Vm + K2 * (q**2 + r**2)/Vm**2
        return D
    
    @staticmethod
    def dynamicEqs(P, D, m, g0, Vm, theta, Sth, Sps, q, r, Cmt, Ne):
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
            Cmt: 
            Ne: # of equations.
        Returns:
            y: Equations for integration.
        """

        # 因变量        
        # y: 
        # index Value
        # 0     vmdot: acceleration of missile 
        # 1     thetadot: first derivative of theta vesus time 
        # 2     psidot: first derivative of psi versus time.
        # 3     pitchdot: pitch rate
        # 4     rdot: yaw rate

        y = np.zeros((5, 1)) # 初始化

        y[0] = (P - D)/m -g0*sin(theta)
        y[1] = (q - cos(theta))/Vm
        y[2] = r/(Vm*cos(theta))
        y[3] = (Sth - q)/Cmt
        y[4] = (Sps - r)/Cmt
        return y

    def odeRK4(self, tEnd, Nt, y0):
        """
        利用四阶龙格库塔法积分运算
        returns:
            y: y[0,:], y[1,:], y[2,:]分别代表Vm、theta、psi
        """
        y = np.zeros((,Nt))
        y[:,0] = y0 # 初值
        t = np.arange(0, tEnd, Nt)
        dt = tEnd/Nt # simulation time interval
        
        for i in range(Nt)
            k1 = dynamicEqs(y[:,i])
            k2 = dynamicEqs(y[:,i] + dt * k1/2)
            k3 = dynamicEqs(y[:,i] + dt * k2/2)
            k4 = dynamicEqs(y[:,i] + dt * k3)
            y[:,i+1] = y[:,i] + dt/6 * (k1 + 2*k2 +2*k3 +k4)            
        return y

    def odeMotions(V, theta, psi, tEnd, Nt):
        """
        对速度方程进行积分得到位置矢量
        Args: 
            V: shape(1, Nt), velocity.
            theta: shape(1, Nt), pitch angle.
            psi: shape(1, Nt), yaw angle.
            tEnd: iteration time.
            Nt: # of time counts.
        returns:
            R: shape(3, Nt), positon vector, [X; Y; Z]
        """        
        dt = tEnd / Nt
        R = np.zeros((3, Nt))
        
        dx = V*cos(theta)*cos(psi)*dt
        dy = V*cos(theta)*sin(psi)*dt
        dz = V*sin(theta)*dt

        R[0, :] += dx
        R[1, :] += dy
        R[2, :] += dz

        R = np.cumsum(R, axis=0)
        return R

    def cmptMass(mMax, mMin, tm, tEnd, Nt):
        """ 
        compute mass, considering mass lost due to oxidation.  
        Args:
            mMax: maximun mass.
            mMin: mass after burn-out.
            tm: work time of motor.
            tEnd: 
            Nt: 
        return:
            m: mass, shape(1, Nt)
        """
        k = (mMax- mMin)/tm
        m = []
        for (i, t) in zip(np.range(Nt), np.arange(0, tEnd, Nt)):
            m[i] = (mMax- k*t) if (t<tm) else mMin
        return np.ndarray(m)

    def cmptDist(Rm, Rt):
        """
        Args:
            Rm: shape(3, Nt), [xm; ym; zm], position vector of missile.
            Rt: shape(3, Nt), [xt; yt; zt], position vector of target.
        return:
            deltaDist: shape(1, Nt), distance of missile and target.
            delta: shape(3, Nt), [deltaX; deltaY; deltaZ], difference of target and missile.
        """
        delta = Rt - Rm

        deltaDist = np.sqrt(np.sum((Rm - Rt)**2, axis=0))
        return deltaDist, delta

    def cmptClsVel(deltaDist, delta, tEnd, Nt):
        """ 
        compute closing velocity
        Args:
            deltaDist: shape(1, Nt), distance of missile and target.
            delta: shape(3, Nt), [deltaX; deltaY; deltaZ], difference of target and missile.
            tEnd:
            Nt:
        return:
            Vc: shape(1, Nt-1), closing velocity.
            CR: shape(3, Nt-1), closing rate, [CRx; CRy; CRz], compute closing rate in X, Y, Z positions.
        """
        dt = tEnd / Nt
        CR = -np.diff(delta, axis=1)/dt

        Vc = np.diff(deltaDist)/dt
        return Vc

    def cmptDesireVel(CR, delta, deltaDist, theta, psi):
        """
        compute desire velocities in position X, Y, Z.
        Args:            
            CR: shape(3, Nt-1), closing rate, [CRx; CRy; CRz], compute closing rate in X, Y, Z positions.
            deltaDist: shape(1, Nt), distance of missile and target.
            delta: shape(3, Nt), [deltaX; deltaY; deltaZ], difference of target and missile.
            theta:
            psi
        return:
            Vdesire: shape(3, Nt-1), desired velocities in the body axis.
            qDesire:
            rDesire:
        """        
        Vdesire = []
        for i in range(Nt-1):
            Vdesire[i] = np.cross(delta[:,i], CR[:,i])/deltaDist[i] 
        Vdesire = np.ndarray(Vdesire)
        qDesire = -np.multiply(sin(psi), Vdesire[0,:]) + np.multiply(cos(psi), Vdesire[1,:])
        rDesire = np.multiply(cos(theta), Vdesire[1,:]) + np.multiply(sin(theta), np.multiply(cos(psi, Vdesire[0,:])+np.multiply(sin(psi,Vdesire[1,:]))))
        return Vdesire, qDesire, rDesire

    def cmptsignals(Cnav, Vc, qDesire, rDesire):
        """ 
        compute the pitch and yaw signals.

        Args:
            Cnav: guidance constants.
            Vc: shape(1, Nt-1), closing velocity.
            qDesire:
            rDesire:
        return:
            Sth: shape(), the pitch(theta) signal.
            Sps: the yaw(psi) signal.
        """
        Sth = Cnav * Vc * qDesire
        Sps = Cnav * Vc * rDesire
        return Sth, Sps


    def simulation(self, parameter_list):
        pass
