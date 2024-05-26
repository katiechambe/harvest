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

        err_str = "No base directory specified. Please specify path to harvest base directory i.e. SetupPaths(/Users/name/Documents/harvest)"
            
        try:
            if basedir == "":
                raise OSError(err_str)
        except OSError:
            raise
        
        self.path_basedir = basedir
        
        self.path_data = self.path_basedir + "/data/"
        self.path_pairs = self.path_data + "pairs/"
        self.path_orbits = self.path_data + "orbits/"
        self.path_misc = self.path_data + "misc/"

    def tng(self, tng_path):
        self.tng_base = tng_path
        self.tng_trees = self.tng_base + "postprocessing/"
        # self.tng_trees = self.tng_base + "postprocessing/"

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
      

