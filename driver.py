import numpy as np
import random


def random_spin(s0, alpha=0.1):
    """Generate a new random spin based on the original one.

    Parameters
    ----------
    s0: np.ndarray

        The original spin to be changed.

    alpha: float

        Larger alpha, larger the change of the spin. Defaults to 0.1.

    Returns
    -------
    np.ndarray

        New updated spin, normalised to 1.

    """
    delta_s = (2 * np.random.random(3) - 1) * alpha
    s1 = s0 + delta_s
    return s1 / np.linalg.norm(s1)


class Driver:
    """Driver class.

    Driver class does not take any input parameters at initialisation. This class
    iterates through random lattice vector changes to minimise the total energy 
    of the system.

    """

    def __init__(self):
        pass

    
    
    def drive(self, system, n, alpha=0.1):

        '''
        system - system
        n - iterations to loop over (potential changen to be mad)s )

        This function  changes one vector at a random lattice point by a random amount
        and compares the change in energy, keeping the lowest energy arrangment and iterating 
        the process n times
        '''

        #obtaining system dimensions
        nx, ny, _ = system.s.array.shape
        
        #counter to break while loop after n interations
        for _ in range(n):
            #determining original system energy
            system_energy = system.energy()

            #choosing a random vector element to change
            rand_i = random.randint(0, nx-1)
            rand_j = random.randint(0, ny-1)
            rand_point = system.s.array[rand_i, rand_j,:].copy()
            rand_s = random_spin(rand_point, alpha=0.1)
            
            #changing random element by a random amount in  copied system 
            system.s.array[rand_i,rand_j,:] = rand_s
            

            #obtaining new system energy
            new_system_energy = system.energy()
    
            
            #comparing system energies and keeping system with lowest energy 
            if new_system_energy < system_energy:
                pass
    
            else:
                system.s.array[rand_i,rand_j,:] = rand_point
        
