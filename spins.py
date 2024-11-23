import numbers
import matplotlib.pyplot as plt
import numpy as np


class Spins:
    """Field of spins on a two-dimensional lattice.

    Each spin is a vector s = (sx, sy, sz). Underlying data structure (``self.array``)
    is a NumPy array (``np.ndarray``) with shape ``(nx, ny, 3)``, where ``nx`` and
    ``ny`` are the number of spins in the x and y directions, respectively, and 3 is for
    three spin (vector) components sx, sy, and sz.

    Parameters
    ----------
    n: Iterable

        Dimensions of a two-dimensional lattice ``n = (nx, ny)``, where ``nx`` and
        ``ny`` is the number of atoms in the x and y directions. Values of ``nx`` and
        ``ny`` must be positive integers.

    value: Iterable

        The value ``(sx, sy, sz)`` initializes all spins in the lattice to have the same
        value. All elements of ``value`` must be real numbers. In general, the values of
        spins do not have the same values, for instance, we can call ``randomise``
        method. Defaults to ``(0, 0, 1)``.

    """

    def __init__(self, n, value=(0, 0, 1)):
        # Checks on input parameters.
        if len(n) != 2:
            raise ValueError(f"Length of n must be 2, not {len(n) = }.")
        if any(i <= 0 or not isinstance(i, int) for i in n):
            raise ValueError("Elements of n must be positive integers.")

        if len(value) != 3:
            raise ValueError(
                f"Length of iterable value must be 3, not {len(value) = }."
            )
        if any(not isinstance(i, numbers.Real) for i in n):
            raise ValueError("Elements of value must be real numbers.")

        self.n = n
        self.array = np.empty((*self.n, 3), dtype=np.float64)
        self.array[..., :] = value  # initialise all spins to have the same value

        # We ensure spins are normalised to 1.
        if not np.isclose(value[0] ** 2 + value[1] ** 2 + value[2] ** 2, 1):
            self.normalise()

    @property
    def mean(self):  
        '''
        Takes no arguments,
        determines mean of each x y and z spin component from every lattice
        point

        '''

        nx = self.n[0]
        ny = self.n[1]

        #summing x y and z components and applying formula
        mean_x = np.sum(self.array[:, :, 0]) / (nx*ny)  
        mean_y = np.sum(self.array[:, :, 1]) / (nx*ny)  
        mean_z = np.sum(self.array[:, :, 2]) / (nx*ny)  

    
        #returning vector of means of components 
        return np.array([mean_x, mean_y, mean_z])


    def __abs__(self):
        '''
        Takes no arguements
        Finds absolute value of each vector at every lattice point
        '''

        A = self.array
        #finding 2x2 lattice dimensions
        i = A.shape[0]
        j = A.shape[1]
        #initializing empty array to append absolute values
        absol = np.empty((i,j))

        #looping over every element and finding absolute value of each vector
        for i1 in range(0,i):

            for j1 in range(0,j):
                #appening to emtpy aray
                absol[i1,j1] = np.linalg.norm(A[i1,j1,:])
        
        return absol[..., np.newaxis]

        
    

    def normalise(self):

        '''
        Takes no arguements 
        Normalizses each vector at every lattice point 
        '''

        self.array = self.array / abs(self)

     

    def randomise(self):

        '''                                                         
        Initialise the lattice with random spins.

        Components of each spin are between -1 and 1: -1 <= si <= 1, and all
        spins are normalised to 1.
        '''
        
                                                                     
        self.array = 2 * np.random.random((*self.n, 3)) - 1
        self.normalise()

    def plot(self):

        '''
        Takes no arguments
        Returns a quiver plot of lattice with spin vectorsd at each lattice
        point

        '''
        #grid size and making grid
        nx, ny, _= self.array.shape
        x = np.arange(nx)
        y = np.arange(ny)
        X, Y = np.meshgrid(x, y)

        #x y and z componenents 
        U = self.array[:, :, 0]  
        V = self.array[:, :, 1] 
        W = self.array[:, :, 2]

        #plotting
        plt.figure()
        plt.quiver(X, Y, U, V,W,scale=None, headwidth=3, headlength=4, headaxislength=4)
        plt.title("2D Lattice of Spins")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.grid(True)
        plt.show()





