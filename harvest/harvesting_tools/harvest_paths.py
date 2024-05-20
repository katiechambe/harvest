""" 
Collection of all paths used in analysis
"""

__author__ = "Katie Chamberlain"
__date__   = "September 2021 - May 2024"

class SetupPaths:

    def __init__(
        self, 
        basedir=""
        ):

        try:
            if basedir == "":
                raise OSError
        except OSError:
            print("Please specify path to harvest base directory."
            print("i.e. SetupPaths(/Users/name/Documents/harvest")
        
        self.path_basedir = basedir
        # self.path_home = self.path_basedir # + "katie/" 
        # self.path_pears = self.path_home + "pears/"

        # # directories for illustris data
        # self.path_illustris = self.path_basedir + "Illustris/"
        # self.path_illustristng = self.path_basedir + "IllustrisTNG/"
        # # group catalog paths
        # self.path_illustrisdark = self.path_illustris + "GroupCatalogsDark/"
        # self.path_illustrishydro = self.path_illustris + "GroupCatalogsHydro/"
        # self.path_tngdark = self.path_illustristng + "TNG100-1-Dark/"
        # self.path_tnghydro = self.path_illustristng + "TNG100-1/"
        # # merger trees
        # self.path_illustrisdark_trees = self.path_illustris + "Illustris-1-Dark-MergerTree/"
        # self.path_illustrishydro_trees = self.path_illustris + "Illustris-1-MergerTree/"
        # self.path_tngdark_trees = self.path_tngdark + "postprocessing/"
        # self.path_tnghydro_trees = self.path_tnghydro + "postprocessing/"
        # # matched catalogs
        # self.path_tngmatch_V = self.path_illustristng + "TNG100-Matched-V/"
        # self.path_tngmatch_N = self.path_illustristng + "TNG100-Matched-Nelson/subhalo_matching_to_dark.hdf5"
      

        # # pears directories
        # self.path_data = self.path_pears + "data/"
        # self.path_groups = self.path_data + "groups/"
        # self.path_subhalos = self.path_data + "subhalos/"
        # self.path_maxmass = self.path_data + "max_masses/"
        # self.path_am_mass = self.path_data + "am_masses/"
        # self.path_snapdata = self.path_data + "snapdata/"
        # self.path_pairs = self.path_data + "pairs/"
        # self.path_median = self.path_data + "median_realization/"
        # self.path_simstars = self.path_data + "simulation_mstar/"
        
        # # pears paper and plots
        # self.path_plots = self.path_pears + "plots/paper1/"
        # self.path_plotdata = self.path_plots + "plotdata/"
        # self.path_paper = self.path_pears + "pears/"