#
#    Vortex Solver ported to python
#
#    author : Vincent Jaunet
#    mail   : vincent.jaunet@ensma.fr
#    date   : december 2019
#
#=====================================================

import numpy as np
import mpmath as mp


class vortex_sheet:
    """
        Class to solve the axisymmetric vortex sheet problem
    """

    def __init__(self,params):
        """ Constructor of the class
            Argument :
                params : a dictionnary with main calculation parameters
                params['Mj'] - jet Mach number
                params['m'] - azimuthal wave number
                params['gama'] - ratio of specific heat
                params['S'] - density ratio  rho_i/rho_o (= To/Ti if ideally expanded)
                params['Rj'] - jet radius

        """
        self.Mj   = params['Mj']        # jet Mach number
        self.m    = params['m']         # azimuthal wave number
        self.gama = params['gama']      # ratio of specific heat
        self.S    = params['S']         # density ratio  rho_i/rho_o (= To/Ti if ideally expanded)
        self.Rj   = params['Rj']        #jet radius

    #======================================================
    def vs_resid(self, k, omega):
        """
        # vs_resid :
        computes the residuals of the Vortex-Sheet
        model for a cylindrical jet

        # Arguments :
        omega  : the pulsation at which we solve the VS
        k      : the wavenumber


         # returns a complex
        """

        u_j = 1                    # jet velocity
        a_j = u_j / self.Mj;       # speed of sound inside jet
        a_o = a_j*np.sqrt(self.S); # speed of sound outside jet

        gama_o = mp.sqrt(k**2 - omega**2)
        gama_i = mp.sqrt(k**2 - 1./self.S*(omega - self.Mj*k)**2)
        if (mp.atan2(gama_i.imag,gama_i.real) < 0):
            gama_i = -gama_i

        Im = mp.besseli(self.m,gama_i*self.Rj);
        Imm = mp.besseli(self.m-1,gama_i*self.Rj);

        Km = mp.besselk(self.m,gama_o*self.Rj);
        Kmm = mp.besselk(int(self.m-1),gama_o*self.Rj);

        resid = 1./(1.-k/omega*self.Mj)**2  + \
                1/self.S*Im/Km*(gama_o*self.Rj*Kmm+self.m*Km)/ \
                (gama_i*self.Rj*Imm - self.m*Im)

        return resid

    #======================================================
    def solve(self, St, k_ini):
        """
            Finding the zeros of the jet VS dispersion relation
            at the frequency St and given k initial guess
        """
        omega = 2*np.pi*St*self.Mj

        k_res = mp.findroot(lambda k: self.vs_resid(k,omega),
                            k_ini, solver='muller')
        return k_res

    #======================================================
    def eigenfunctions(self,St,k):
        r_o=np.linspace(0.5,10.0,200)
        r_i=np.linspace(0,0.5,50)

        omega=2.0*np.pi*St*self.Mj;

        gamai=mp.sqrt(k**2-(omega-self.Mj*k)**2);
        gamao=mp.sqrt(k**2-omega**2)

        C=mp.besselk(0,gamao*0.5)/\
            mp.besseli(0,gamai*0.5)

        besseli_array = np.frompyfunc(mp.besseli, 2, 1)
        besselk_array = np.frompyfunc(mp.besselk, 2, 1)

        p_i = C*besseli_array(0,gamai*r_i)
        p_o =   besselk_array(0,gamao*r_o)
        p_i=np.array(p_i.tolist(),dtype=np.complex64)
        p_o=np.array(p_o.tolist(),dtype=np.complex64)

        # streamwise velocity
        #u_i = -k/(M*k-omega)*p_i;
        #u_o = k/omega*p_o;

        # radial velocity
        #v_i = 1i/(M*k-omega)*diff(p_i)./diff(r_i);
        #v_o = -1i/omega*diff(p_o)./diff(r_o);

        return r_i, p_i, r_o, p_o

        #======================================================
    if __name__ == "__main__":

        import vs_solver
        param={}
        param['Mj']  =1.2
        param['gama']=1.4
        param['m'] = 0
        param['S'] = 1
        param['Rj'] = 0.5
        vs=vs_solver.vortex_sheet(param)

        St = 0.4
        k_ini = 2*np.pi*St*1.2/param['Mj']-1.0j
        k_kh = vs.solve(St,k_ini)
        print(k_kh)
        # expected value k_kh = 2.7317 - 1.9626i
