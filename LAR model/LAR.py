from math import *
import numpy as np

class larCalc(object):
    def __init__(self, parameter_list):
        # m = 500    # kg         
        # P = 0.0    # N
        # g0 = 9.8066         # Acceleration due to gravity [m/s^2]

        # # SET INTEGRATION PARAMETERS 
        # tEnd = 20000        #timeout value

    @staticmethod
    def calcAlphaBeta(P, thetaDot, psiDot, *args):
        """
        Args:
            P: Thrust.
            thetaDot: first derivative of theta vesus time.
            psiDot: first derivative of psi versus time.
            m: mass of missile.
            g0: Gravity acceleration.
            Vm: velocity of missile
            theta: 弹道倾角
            Q: 动压
            S: reference area
            Cy: 
            Cz:             
        returns:
            alpha: attack angle
            beta: slideslip angle
        """
        m, g0, Vm, theta, Q, S, Cy, Cz = args
        alpha = (m * Vm * thetaDot + m * g0 * cos(theta)) / (P + 57.3*Q*S*Cy)
        beta = (m * Vm * psiDot *cos(theta)) / (P - 57.3*Q*S*Cz)
        return alpha, beta

    @staticmethod
    def calcAerodynmcs():
        """
        计算气动力

        Args:
            Q: 动压
            S: reference area
            Cx0:
            Cx:
            Cy: 
            Cz: 
            alphaSum:
        Returns:            
            Fx, Fy, Fz: 空气动力在速度坐标系坐标轴上的投影
        """
        Fx = Q*S*(Cx0 + Cx*(alphaSum*57.3)**2) 
        Fy = 57.3*Q*S*Cy
        Fz = 57.3*Q*S*Cz
            
    
    @staticmethod
    def dynamicEqs(y, params):
        """
        导弹的动力学方程组
        
        Returns:
            y1

        """

        # calculate alpha, beta
        alpha, beta = calcAlphaBeta(P, thetaDot, psiDot, params)

        # calculate aerodynamics force
        Fx, Fy, Fz = calcAerodynmcs()

        # 因变量
        Vm = y[0]
        theta = y[1]
        psi = y[2]
        y1 = np.zeros((3, 1))
        y1[0] = P/m * cos(alpha) * cos(beta) - Fx/m -g0*sin(theta)
        y1[1] = P/(m*Vm)*sin(alpha) - Fy/(m * Vm) - g/Vm*cos(theta)
        y1[2] = -P/(m*Vm*cos(theta))*cos(alpha)*sin(beta) + Fz/(m*Vm*cos(theta))
        return y1

    def odeRK4(self, tEnd, Nt, y0):
        """
        利用四阶龙格库塔法积分运算
        returns:
            y: y[0,:], y[1,:], y[2,:]分别代表Vm、theta、psi
        """
        y = np.zeros((3,Nt))
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
            V: shape(1, Nt)
            theta: shape(1, Nt)
            psi: shape(1, Nt)
            tEnd: iteration time
            Nt: # of time 
        returns:
            R: 位置矢量, shape(3, Nt), [X; Y; Z]
        """
        
        dt = tEnd / Nt
        R = np.zeros((3, Nt))
        
        dx = V*cos(theta)*cos(psi)*dt
        dy = V*sin(theta)*dt
        dz = V*cos(theta)*sin(psi)*dt

        R[0, :] += dx
        R[1, :] += dy
        R[2, :] += dz

        R = np.cumsum(R, axis=0)
        return R

    def computeDist(Rm, Rt):
        """
        Args:
            Rm: [xm, ym, zm], position of missile.
            Rt: [xt, yt, zt], position of target.
        return:
            deltaDist: distance of missile and target.
        """
        deltaDist = np.sqrt(np.sum((Rm - Rt)**2, axis=0))
        return deltaDist

    def computeClsVel(deltaDist):
        """ 
        compute closing velocity
        Args:
            deltaDist: distance of missile and target.
        return:
            Vc: shape(1, Nt-1), closing velocity.
        """
        Vc = np.diff(deltaDist)
        return Vc

    def simulation(self, parameter_list):
        pass
