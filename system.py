import numpy as np
import matplotlib.pyplot as plt



class System:
    """System object with the spin configuration and simulation parameters.

    Parameters
    ----------
    s: mcsim.Spins

        Two-dimensional spin field.

    B: Iterable(float)

        External magnetic field (vector), length 3.

    K: numbers.Real

        Uniaxial anisotropy constant.

    u: Iterable(float)

        Uniaxial anisotropy axis, length 3. If ``u`` is not normalised to 1, it
        will be normalised before the calculation of uniaxial anisotropy energy.

    J: numbers.Real

        Exchange energy constant.

    D: numbers.Real

        Dzyaloshinskii-Moriya energy constant.

    """

    def __init__(self, s, B, K, u, J, D):
        self.s = s
        self.J = J
        self.D = D
        self.B = B
        self.K = K
        self.u = u

    def energy(self):
        """Total energy of the system.

        The total energy of the system is the sum of all individual energy terms.

        Returns
        -------
        float

            Total energy of the system.

        """
        return self.zeeman() + self.anisotropy() + self.exchange() + self.dmi()

    def zeeman(self):

        '''
        Takes no arguements
        Calculates the zeeman energy of the vector lattice
        
        '''
        #sum of each spin component x, y, z at each point
        sum_x = np.sum(self.s.array[:, :, 0])  
        sum_y = np.sum(self.s.array[:, :, 1])  
        sum_z = np.sum(self.s.array[:, :, 2]) 
        
        #array of sums 
        sum_spins = np.array([sum_x, sum_y, sum_z])

        #zeeman formula
        Z = - np.dot(sum_spins,self.B)

        return Z

    def anisotropy(self):
        """
        Takes no arguements
        Calculates anisotropy energy of the vector lattice


        """
        #normalizing vector to manitude 1
        u_mag = np.linalg.norm(self.u)
        u_norm = self.u/u_mag

        s = 0
        for i in range(self.s.n[0]):
            for j in range(self.s.n[1]):
                s += np.dot(self.s.array[i, j, :], u_norm)**2

        return -self.K * s
       


    def exchange(self):
        """
        Takes no arguements
        calculates the exchange energy of the vector lattice 
        
        """
        #finding loop limits
        nx, ny, _ = self.s.array.shape
        
        #initializing value allowing sum over loops
        ex = 0.0

        #looping over each element and dotting with neighbouring and summing 
        for i in range(nx):
            for j in range(ny - 1):
                ex += np.dot(self.s.array[i, j,:], self.s.array[i, j + 1,:])
    
        for j in range(ny):
            for i in range(nx - 1):
                ex += np.dot(self.s.array[i, j,:], self.s.array[i + 1, j,:])

        #J constant 
        ex = -self.J*ex

        return ex

    def dmi(self):
        '''
        Takes no arguements
        Calculates the DMI energy of the vector lattice
        
        '''
        
        #finding loop limits
        nx, ny, _ = self.s.array.shape
        #initializing value allowing sum over loops
        dmi = 0

        #looping over each element and crossing with neighbouring and summing 
        for i in range(nx):
            for j in range(ny - 1):
                dmi += np.dot([self.D,0,0],(np.cross(self.s.array[i, j,:], self.s.array[i, j + 1,:])))

        #D vector pointing towards neighbour 
        
        for j in range(ny):
            for i in range(nx - 1):
                dmi += np.dot([0,self.D,0],(np.cross(self.s.array[i, j,:], self.s.array[i + 1, j,:])))



        return dmi

    

