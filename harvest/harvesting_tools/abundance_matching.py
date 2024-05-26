""" 
Creates abundance matched realizations for dark matter halos

Usage:
------
To get stellar mass:
    AbundanceMatching(maxmass,redshift,#samples).stellar_mass(med)
    -- NOTE:
        - med=False gives a vector of stellar masses with some spread and
            length samples
        - med=True yields one individual value for the given max mass and
            redshift and is calculated by setting errors to 0

Details:
--------
Follows Moster, Naab, and White (2012)
https://arxiv.org/pdf/1205.5807.pdf
equations 2,11-14, and table 1
"""

__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "May 2019 - edited Oct. 2021"

from numpy.random import normal
from numpy import where

class AbundanceMatching:
    def __init__(self, maxmass, redshift, samples):
        """
        Samples from gaussian distributed abundance matching relation 

        Parameters:
        -----------
        maxmass: float
            the dark matter maximum halo mass ever achieved
            (in Msun, NOT 1e10Msun)
        redshift:
            redshift at the present snapshot
            NOTE - not the snapshot at max mass
        samples: int
            number of realizations to sample from relationship  
        """
        
        self.maxmass = maxmass
        self.z = redshift
        self.samples = samples

        # values from table 1
        self.M10vec = [11.59, 0.236]
        self.M11vec = [1.195, 0.353]
        self.N10vec = [0.0351, 0.0058]
        self.N11vec = [-0.0247, 0.0069]
        self.beta10vec = [1.376, 0.153]
        self.beta11vec = [-0.826, 0.225]
        self.gamma10vec = [0.608, 0.059]
        self.gamma11vec = [0.329, 0.173]

    def getvals(self, vec, med=False):
        """ 
        if asking for median (med=True):
            then the median value is returned
        if med=False, 
            return a vector of length self.samples of 
            values with gaussian error
        """
        if med:
            return normal(vec[0], 0)
        else:
            return normal(vec[0], vec[1], self.samples)

    def func(self, x, dx):
        """of form x + dx*(z/z+1)"""
        return x + dx*(self.z/(1+self.z))

    def mass_ratio(self, med=False):
        """
        stellar to halo mass ratio

        Returns:
        --------
            Stellar mass to halo mass ratio

        eq.2 in Moster
        """
        M10 = self.getvals(self.M10vec, med)
        M11 = self.getvals(self.M11vec, med)
        N10 = self.getvals(self.N10vec, med)
        N11 = self.getvals(self.N11vec, med)
        beta10 = self.getvals(self.beta10vec, med)
        beta11 = self.getvals(self.beta11vec, med)
        gamma10 = self.getvals(self.gamma10vec, med)
        gamma11 = self.getvals(self.gamma11vec, med)
    
        logM = self.func(M10, M11)
        Nany = self.func(N10, N11)
        Nmed = self.func( self.getvals(self.N10vec, "True"), self.getvals(self.N11vec, "True") )
        Npos = where(Nany < 0, Nmed, Nany)
        beta = self.func(beta10, beta11)
        gamma = self.func(gamma10, gamma11)

        M1 = 10**logM
        A = (self.maxmass/M1)**(-beta)
        B = (self.maxmass/M1)**(gamma)
        SHMratio = 2*Npos*(A+B)**-1
        return SHMratio

    def stellar_mass(self,med=False):
        """ 
        stellar mass

        Returns:
        --------
            stellar mass in Msun
        """
        return self.maxmass*self.mass_ratio(med)

 
