####!/usr/bin/env python
""" routines for reading subfind data from cosmo sims.
    
    Example Usage:
    
    Dependencies:
      hdf5lib.py


    Notice:
    -------
    Katie Chamberlain modified this code slightly in 2021 to resolve 
    deprecation
"""

__author__ = "Mark Vogelsberger, Paul Torrey and contributing authors"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Mark Vogelsberger, Paul Torrey and contributing authors"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@mit.harvard.edu"
__status__ = "Beta -- forever."


import os
import sys
import numpy as np
import hdf5libPy3 as hdf5lib
import namingPy3 as naming


####################
#SUBHALO DATABLOCKS#
####################
#descriptions of subhalo datablocks -> add new datablocks here!
sub_datablocks = {"SubhaloLen":["INT",1],
                  "SubhaloMass":["FLOAT",1],
                  "SubhaloMassinRad":["FLOAT",1],
                  "SubhaloPos":["FLOAT",3],
                  "SubhaloVel":["FLOAT",3],
                  "SubhaloLenType":["INT",6],
                  "SubhaloMassType":["FLOAT",6],
                  "SubhaloCM":["FLOAT",3],
                  "SubhaloSpin":["FLOAT",3],
                  "SubhaloVelDisp":["FLOAT",1],
                  "SubhaloVmax":["FLOAT",1],
                  "SubhaloVmaxRad":["FLOAT",1],
                  "SubhaloHalfmassRad":["FLOAT",1],
                  "SubhaloHalfmassRadType":["FLOAT",6],
                  "SubhaloMassInRadType":["FLOAT", 6],
                  "SubhaloMassInRad":["FLOAT",1],
                  "SubhaloMassInHalfRadType":["FLOAT", 6],
                  "SubhaloMassInHalfRad":["FLOAT", 1],
                  "SubhaloIDMostbound":["ID",1],
                  "SubhaloGrNr":["INT",1],
                  "SubhaloParent":["INT",1],
                  "SubhaloSFR":["FLOAT",1],
                  "SubhaloSFRinRad":["FLOAT",1],
                  "SubhaloGasMetallicity":["FLOAT",1],
                  "SubhaloGasMetallicitySfr":["FLOAT",1],
                  "SubhaloStarMetallicity":["FLOAT",1],
                  "SubhaloGasMetalFractions":["FLOAT",9],
                  "SubhaloGasMetalFractionsSfr":["FLOAT",9],
                  "SubhaloGasMetalFractionsSfrWeighted":["FLOAT",9],
                  "SubhaloStarMetalFractions":["FLOAT",9],
                  "SubhaloStarMetallicityHalfRad":["FLOAT",1],
                  "SubhaloBHMass":["FLOAT",1],
                  "SubhaloBHMdot":["FLOAT",1],
                  "SubhaloStellarPhotometricsMassInRad":["FLOAT",1],
                  "SubhaloStellarPhotometrics":["FLOAT",8]}  #band luminosities: U, B, V, K, g, r, i, z

##################
#GROUP DATABLOCKS#
##################
#descriptions of subhalo datablocks -> add new datablocks here!
#format -> "HDF5_NAME":["DATATYPE", DIMENSION]
grp_datablocks = {"GroupLen":["INT",1],
                  "GroupMass":["FLOAT",1],
                  "GroupPos":["FLOAT",3],
                  "GroupVel":["FLOAT",3],
                  "GroupLenType":["INT",6],
                  "GroupMassType":["FLOAT",6],
                  "Group_M_Mean200":["FLOAT",1],
                  "Group_R_Mean200":["FLOAT",1],
                  "Group_M_Crit200":["FLOAT",1],
                  "Group_R_Crit200":["FLOAT",1],
                  "Group_M_TopHat200":["FLOAT",1],
                  "Group_R_TopHat200":["FLOAT",1],
                  "Group_M_Crit500":["FLOAT",1],
                  "Group_R_Crit500":["FLOAT",1],
                  "GroupNsubs":["INT",1],
                  "GroupFirstSub":["INT",1],
                  "GroupSFR":["FLOAT",1],
                  "GroupGasMetallicity":["FLOAT",1],
                  "GroupStarMetallicity":["FLOAT",1],
                  "GroupGasMetalFractions":["FLOAT",9],
                  "GroupStarMetalFractions":["FLOAT",9],
                  "GroupBHMass":["FLOAT",1],
                  "GroupBHMdot":["FLOAT",1],
                  "GroupFuzzOffsetType":["INT64",6]}

galprop_datablocks = {
                        "nh_mass_inmidrad":         ["FLOAT",1, "SubfindNHIMassInThreeRad"],
                        "nh_mass_inrad":            ["FLOAT",1, "SubfindNHIMassInTwoRad"],
                        "nh_mass_intworad":         ["FLOAT",1, "SubfindNHIMassInFourRad"],
                        "nh_sfrmass":               ["FLOAT",1, "SubfindNHISFMass"],
                        "nh_totmass":               ["FLOAT",1, "SubfindNHITotMass"],
                        "stellar_mag_inrad":        ["FLOAT",8, "SubfindStellarMagInTwoRad"],
                        "stellar_mag_intworad":     ["FLOAT",8, "SubfindStellarMagInFourRad"],
                        "stellar_metallicity_inrad":["FLOAT",1, "SubfindStellarMetInTwoRad"],
                        "wind_totmass":             ["FLOAT",1, "SubfindTotWindMass"]
                    }




class subfind_catalog:
    def __init__(self, basedir, snapnum, long_ids=False, double_output=False, grpcat=True, subcat=True, name="fof_subhalo_tab", keysel=[]):

        if long_ids: self.id_type = np.uint64
        else: self.id_type = np.uint32
        if double_output: self.double_type = np.float32
        else: self.double_type = np.float64

        filenum = 0
        doneflag = False
        skip_gr = 0
        skip_sub = 0
        vardict = {}
        if keysel is None:
            keysel = grp_datablocks.items()

        while not doneflag:
            self.filebase, curfile = naming.return_subfind_filebase(basedir, snapnum, name, filenum)
            self.firstfile = curfile

            f=hdf5lib.OpenFile(curfile)
            ngroups = hdf5lib.GetAttr(f, "Header", "Ngroups_ThisFile")
            nsubs = hdf5lib.GetAttr(f, "Header", "Nsubgroups_ThisFile")
            nfiles = hdf5lib.GetAttr(f, "Header", "NumFiles")
            if filenum == 0:
                self.ngroups = hdf5lib.GetAttr(f, "Header", "Ngroups_Total")
                self.nids = hdf5lib.GetAttr(f, "Header", "Nids_Total")
                self.nsubs = hdf5lib.GetAttr(f, "Header", "Nsubgroups_Total")
                self.redshift = hdf5lib.GetAttr(f, "Header", "Redshift")
                self.boxsize = hdf5lib.GetAttr(f, "Header", "BoxSize")
                #GROUPS
                if grpcat:
                    for key in keysel:
                        if hdf5lib.Contains(f, "Group", key):
                            val = grp_datablocks[key]
                            type = val[0]
                            dim = val[1]
                            
                            # Notice: added if and elif to get rid of deprecation error:
                            # Previously: 
                            # if (type=='FLOAT'):
                            #    vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((self.double_type,dim)))
                            if (type=='FLOAT') and (dim ==1):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((self.double_type)))
                            elif (type=='FLOAT'):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((self.double_type,dim)))

                            if (type=='INT') and (dim==1):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((np.int32)))
                            elif (type=='INT'):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((np.int32,dim)))

                            if (type=='INT64') and (dim==1):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((np.int64)))
                            elif (type=='INT64'):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((np.int64,dim)))

                            if (type=='ID') and (dim==1):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((self.id_type)))
                            if (type=='ID'):
                                vars(self)[key]=np.empty(self.ngroups, dtype=np.dtype((self.id_type,dim)))
                            vardict[key]=vars(self)[key]


                #SUBHALOS
                if subcat:
                    for key in keysel:
                        if hdf5lib.Contains(f, "Subhalo", key):
                            val = sub_datablocks[key]
                            type = val[0]
                            dim = val[1]
                            # Notice: added if and elif to get rid of deprecation error
                            if (type=='FLOAT') and (dim ==1):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((self.double_type)))
                            elif (type=='FLOAT'):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((self.double_type,dim)))

                            if (type=='INT') and (dim==1):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((np.int32)))
                            elif (type=='INT'):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((np.int32,dim)))

                            if (type=='INT64') and (dim==1):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((np.int64)))
                            elif (type=='INT64'):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((np.int64,dim)))

                            if (type=='ID') and (dim==1):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((self.id_type)))
                            if (type=='ID'):
                                vars(self)[key]=np.empty(self.nsubs, dtype=np.dtype((self.id_type,dim)))
                                
                            vardict[key]=vars(self)[key]

            #GROUPS
            if grpcat:
                if ngroups > 0:
                    for key in keysel:
                        if hdf5lib.Contains(f, "Group", key):
                            val = grp_datablocks[key]
                            type = val[0]
                            dim = val[1]
                            a=hdf5lib.GetData(f, "Group/"+key)
                            if dim==1:
                                vardict[key][skip_gr:skip_gr + ngroups]=a[:]
                            else:
                                for d in range(0,dim):
                                    vardict[key][skip_gr:skip_gr + ngroups,d]=a[:,d]

                    skip_gr += ngroups
            #SUBHALOS
            if subcat:
                if nsubs > 0:
                    for key in keysel:
                        if hdf5lib.Contains(f, "Subhalo", key):
                            val = sub_datablocks[key]
                            type = val[0]
                            dim = val[1]
                            a=hdf5lib.GetData(f, "Subhalo/"+key)
                            if dim==1:
                                vardict[key][skip_sub:skip_sub + nsubs]=a[:]
                            else:
                                for d in range(0,dim):
                                    vardict[key][skip_sub:skip_sub + nsubs,d]=a[:,d]

                    skip_sub += nsubs

            f.close()

            filenum += 1
            if filenum == nfiles: doneflag = True



    def list_blocks(self, parttype=-1, verbose=False):
        curfile = self.firstfile
        if not os.path.exists(curfile):
            print("file not found:", curfile)
            sys.exit()

        f=hdf5lib.OpenFile(curfile)
        iter = it=sub_datablocks.__iter__()
        next = iter.next()
        while (1):
            print(next)
            if (hdf5lib.Contains(f,"Subhalo",next)):
                print("Subhalo: "+next)
                sys.stdout.flush()
            try:
                next=iter.next()
            except StopIteration:
                break
        f.close()



class galprop_catalog:
    def __init__(self, basedir, snapnum, keysel=None, long_ids=False):
        
        if long_ids: id_type = np.uint64
        else: id_type = np.uint32
        
        vardict = {}
        if keysel is None:
            keysel = galprop_datablocks.items()
        
        file=naming.return_galprop_file(basedir, snapnum)
        if os.path.exists(file):
            f=hdf5lib.OpenFile(file, mode='r')
            for key in keysel:
                if hdf5lib.Contains(f, "", key):
                    val = galprop_datablocks[key]
                    type = val[0]
                    dim = val[1]
                    vars(self)[key]=np.array(hdf5lib.GetData(f, key)[:])
            f.close()
        else:
            print("Galprop File Not Found...")

        
        
        








# TODO: Improve this.
def get_offsets(cat, part_types=[0, 1, 4, 5], snap=None, run=None):
    if snap and run:
        group_file = "/n/ghernquist/Illustris/Runs/%s/postprocessing/offsets/snap_offsets_group_%s.hdf5" % (run, snap)
        halo_file = "/n/ghernquist/Illustris/Runs/%s/postprocessing/offsets/snap_offsets_subhalo_%s.hdf5" % (run, snap)
        if os.path.isfile(group_file) and os.path.isfile(halo_file):
            print("READSUBF: found pretabulated offsets to read")
            f = hdf5lib.OpenFile(group_file)
            group_offsets = hdf5lib.GetData(f, "Offsets")[:]
            f.close()

            f = hdf5lib.OpenFile(halo_file)
            halo_offsets = hdf5lib.GetData(f, "Offsets")[:]
            f.close()

            return np.array(group_offsets), np.array(halo_offsets)

        else:
          # /n/hernquistfs3/IllustrisTNG/Runs/L75n910TNG/postprocessing/offsets/
            group_file = "/n/hernquistfs3/IllustrisTNG/Runs/%s/postprocessing/offsets/offsets_%s.hdf5" % (run, str(snap).zfill(3)) 	
#            sys.exit()
            if os.path.isfile(group_file):
                f = hdf5lib.OpenFile(group_file)
                group_offsets = np.copy(hdf5lib.GetData(f, "Group/SnapByType"))
                halo_offsets  = np.copy(hdf5lib.GetData(f, "Subhalo/SnapByType"))
                return group_offsets, halo_offsets
     
   
    GroupOffset = np.zeros((cat.ngroups, 6), dtype="int64")
    HaloOffset  = np.zeros((cat.nsubs, 6), dtype="int64")

    for parttype in part_types:
        print("Calculating offsets for PartType: %d" % parttype)
        k = 0
        for i in range(0, cat.ngroups):
                    if i > 0:
                           GroupOffset[i, parttype] = GroupOffset[i-1, parttype] + cat.GroupLenType[i-1, parttype]
                    if cat.GroupNsubs[i] > 0:
                            HaloOffset[k, parttype] = GroupOffset[i, parttype]
                            k += 1
                            for j in range(1, cat.GroupNsubs[i]):
                                    HaloOffset[k, parttype] =  HaloOffset[k-1, parttype] + cat.SubhaloLenType[k-1, parttype]
                                    k += 1
    if k != cat.nsubs:
        print("READHALO: problem with offset table", k, cat.nsubs)
        sys.exit()

    return np.array(GroupOffset), np.array(HaloOffset)


def subhalo_offsets(snap = 135, run='Illustris-1'):
    snaptag=str(snap)
    f=hdf5lib.OpenFile('/n/ghernquist/Illustris/Runs/'+run+'/postprocessing/offsets/snap_offsets_subhalo_'+snaptag+'.hdf5', mode ='r' )
    data=hdf5lib.GetData(f, "Offsets")[:]
    f.close()
    return np.array(data)

def subhalo_insitu_fraction(snap = 135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile('/n/ghernquist/Illustris/Runs/Illustris-1/postprocessing/InSituFraction/insitu_stellar_fraction_'+snaptag+'.hdf5', mode ='r' )
    data=hdf5lib.GetData(f, "InSitu")[:]
    f.close()
    return np.array(data)

def subhalo_stellar_vel_disp(base, snap = 135, which="StellarVelDisp_HalfMassRad"):
    #snaptag='000'+str(snap)
    #snaptag=snaptag[-3:]
    snaptag=str(snap)
    print(base+'/postprocessing/stellar_vel_disp/stellar_vel_disp_'+snaptag+'.hdf5')
    f=hdf5lib.OpenFile(base+'/postprocessing/stellar_vel_disp/stellar_vel_disp_'+snaptag+'.hdf5', mode ='r' )
    delta=hdf5lib.GetData(f, "Subhalo/"+which)[:]
    f.close()
    return np.array(delta)

def subhalo_gas_z_grad(base, snap = 135, which="GradMetallicity_5"):
    snaptag=str(snap)
    file = base+'/postprocessing/gas_metallicity/gas_metallicity_info_'+snaptag+'.hdf5'
    print(file)
    f=hdf5lib.OpenFile(file, mode ='r' )
    data=hdf5lib.GetData(f, "Subhalo/"+which)[:]
    f.close()
    return np.array(data)

def subhalo_gas_kinematics(base, snap = 135, which="v_5"):
    snaptag=str(snap)
    file = base+'/postprocessing/gas_kinematics/gas_kinematic_info_'+snaptag+'.hdf5'
    print(file)
    f=hdf5lib.OpenFile(file, mode ='r' )
    data=hdf5lib.GetData(f, "Subhalo/"+which)[:]
    f.close()
    return np.array(data)

def subhalo_overdensity(base, snap = 135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/environment/environment_'+snaptag+'.hdf5', mode ='r' )
    delta=hdf5lib.GetData(f, "delta")[:]
    f.close()
    return np.array(delta)

def insitu_mass(base, snap=135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/StellarAssembly/galaxies_'+snaptag+'.hdf5', mode ='r' )
    delta=hdf5lib.GetData(f, "StellarMassInSitu")[:]
    f.close()
    return np.array(delta)

def exsitu_mass(base, snap=135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/StellarAssembly/galaxies_'+snaptag+'.hdf5', mode ='r' )
    delta=hdf5lib.GetData(f, "StellarMassExSitu")[:]
    f.close()
    return np.array(delta)

def mass_from_mergers(base, snap=135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/StellarAssembly/galaxies_'+snaptag+'.hdf5', mode ='r' )
    delta=hdf5lib.GetData(f, "StellarMassFromMergers")[:]
    f.close()
    return np.array(delta)


def mass_from_major_mergers(base, snap=135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/MergerHistory/merger_history_'+snaptag+'.hdf5', mode ='r' )
    data=hdf5lib.GetData(f, "StellarMassFromMajorMergers")[:]
    f.close()
    return np.array(data)

def mass_from_minor_mergers(base, snap=135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/MergerHistory/merger_history_'+snaptag+'.hdf5', mode ='r' )
    data=hdf5lib.GetData(f, "StellarMassFromMinorMergers")[:]
    f.close()
    return np.array(data)

def number_of_major_mergers(base, snap=135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/MergerHistory/merger_history_'+snaptag+'.hdf5', mode ='r' )
    data=hdf5lib.GetData(f, "NumMajorMergersTotal")[:]
    f.close()
    return np.array(data)


def number_of_minor_mergers(base, snap=135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    f=hdf5lib.OpenFile(base+'/postprocessing/MergerHistory/merger_history_'+snaptag+'.hdf5', mode ='r' )
    data=hdf5lib.GetData(f, "NumMinorMergersTotal")[:]
    f.close()
    return np.array(data)



def subhalo_circularities(base, snap = 135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    print(base+'/postprocessing/circularities/circularities_'+snaptag+'.hdf5')
    f=hdf5lib.OpenFile(base+'/postprocessing/circularities/circularities_'+snaptag+'.hdf5', mode ='r' )
    data=np.array(hdf5lib.GetData(f, "CircAbove05Frac")[:])
    data=np.reshape(data, -1)
    f.close()
    return data



def subhalo_stellar_metallicities(snap = 135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    file='/n/ghernquist/Illustris/Runs/Illustris-1/postprocessing/galprop/galprop_'+snaptag+'.hdf5'
    if os.path.exists(file):
        f=hdf5lib.OpenFile('/n/ghernquist/Illustris/Runs/Illustris-1/postprocessing/galprop/galprop_'+snaptag+'.hdf5', mode='r')
        data=np.array(hdf5lib.GetData(f, "stellar_metallicity_inrad")[:])
        f.close()
    else:
        data = None
    return data

def subhalo_stellar_age(snap = 135):
    snaptag='000'+str(snap)
    snaptag=snaptag[-3:]
    file='/n/ghernquist/Illustris/Runs/Illustris-1/postprocessing/galprop/galprop_'+snaptag+'.hdf5'
    if os.path.exists(file):
        f=hdf5lib.OpenFile('/n/ghernquist/Illustris/Runs/Illustris-1/postprocessing/galprop/galprop_'+snaptag+'.hdf5', mode='r')
        data=np.array(hdf5lib.GetData(f, "stellar_age_inrad")[:])
        f.close()
    else:
        data = None
    return data


def subhalo_petrosian_radius(snap = 135):
    file = '/n/ghernquist/Illustris/Runs/Illustris-1/postprocessing/PhotometricMorphologies/nonparmorphs_iSDSS_135.hdf5'
    if os.path.exists(file):
        f=hdf5lib.OpenFile(file,mode='r')
        data0=np.array(hdf5lib.GetData(f, "RP_cam0")[:])
        data = np.zeros( (4 , data0.shape[0]) )

        data[1,:] = np.array(hdf5lib.GetData(f, "RP_cam1")[:])
        data[2,:] = np.array(hdf5lib.GetData(f, "RP_cam2")[:])
        data[3,:] = np.array(hdf5lib.GetData(f, "RP_cam3")[:])

        data = np.median( data, axis=0 )

        f.close()
    else:
        data = None
    return data


    

