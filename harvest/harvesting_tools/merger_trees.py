import numpy as np
from harvesting_tools.readtreeHDF5_public import TreeDB


class TraceMergerTree:

    def __init__(
        self, 
        treepath,
        snapshot,
        subfindID
        ):
        """
        Identifies and pulls merger tree for a single subhalo

        Parameters
        ----------
        snapshot: int
            the number of the snapshot with the corresponding subhalo ID
        subfindID: int
            the ID number of the subhalo at the corresponding snapshot
        """
        self.snapshot = snapshot
        self.subfindID = subfindID

        self.treeDirectory = treepath

        tree = TreeDB(self.treeDirectory)
        pastbranch = tree.get_main_branch( 
            self.snapshot, 
            self.subfindID
            # keysel=['SnapNum', 'SubhaloMass', 'SubhaloPos', 'SubhaloVel', 'SubhaloID', 'SubfindID']
            )
        futurebranch = tree.get_future_branch(
                               self.snapshot,
                               self.subfindID)

        self.pastbranch = pastbranch
        self.futurebranch = futurebranch
        
        self.pastkeys = np.array(list(self.pastbranch.__dict__.keys()))
        self.futurekeys = np.array(list(self.futurebranch.__dict__.keys()))
        
        self.fullbranch = {}
        for key in self.pastkeys[np.isin(self.pastkeys,self.futurekeys)]:
        #print(type(tree1.futurebranch.__getattribute__(key)))
                self.fullbranch[key] = np.concatenate([self.futurebranch.__getattribute__(key)[:-1],self.pastbranch.__getattribute__(key)])[::-1]
        
    @property
    def maxmass(self):
        """
        Max mass of the subhalo
        -- note: this only considers current and previous snapshots --

        Parameters:
        -----------
        None

        Outputs:
        --------
        maxmass: float
            the maximum mass previously achieved by a subhalo
        maxsnap: int
            the snapshot at which max mass occurs 
        maxredshift: float
            the maximum redshift at which max mass occurs
        """
        maxmass = max(self.masses_phys)
        maxmass_mask =  max(self.masses_phys)==self.masses_phys
        maxsnap = self.snaps[maxmass_mask][0]
        return maxmass, maxsnap


