import os
import sys
import h5py
import numpy as np

from harvesting_tools.harvest_paths import SetupPaths
import illustris_tools.readsubfHDF5_Py3 as readSub

from utils_orbs.merger_trees import TraceMergerTree
from utils_orbs.readMergerHistoryFiles import FindHistory
from utils_orbs.vectorCorrection import vectorCorrection as vector


path_to_harvest = "/xdisk/gbesla/katiechambe/harvest"
paths = SetupPaths(path_to_harvest)

snapshot = int(sys.argv[1])
masstype = str(sys.argv[2])
pairtype = str(sys.argv[3])

print("Read in all packages and defined snapshot")

full_snaps = np.arange(0,100,1)
little_h = 0.6774

# check if file already exists.
# if file exists, finish running
print(paths.path_data)
if os.path.exists(f"{paths.path_data}orbits/{masstype}mass_{pairtype}_{snapshot}.hdf5"):
    print ("file already exists")

# if file doesn't exist, create it
else:
    # ----------------------------------------------------
    # import pair data
    f = h5py.File(f"{paths.path_data}pairs/{masstype}mass_{pairtype}_{snapshot}.hdf5", 'r')
    pairs = {}
    for key, val in f.items():
        if key != "Header":
            pairs[key] = np.array(val) 
    f.close()
    numpairs = len(pairs['Sub1 ID'])

    # ----------------------------------------------------
    # import snapshot data to convert comoving units
    f = h5py.File(f"{paths.path_data}misc/snapshot_data.hdf5","r")
    convert = {}
    for key,val in f.items():
        convert[key] = np.array(val)
    f.close()
    
    # heading info:
    redshift_snap = convert['Redshift'][snapshot]
    scale_snap = convert['Scale'][snapshot]

    # ----------------------------------------------------
    # initialize data structures
    
    # 1D (per pair) - 
    groupnum  = np.zeros( numpairs, dtype="int32") # group number of primary
    id1, id2  = np.zeros((2, numpairs), dtype="int32") # subhalo ID 1 & 2
    mass1, mass2 = np.zeros((2, numpairs)) # subhalo mass 1 & 2
    stell1, stell2 = np.zeros((2, numpairs)) # stellar mass 1 & 2
    smratio  = np.zeros(numpairs) # stellar mass ratio (mass2/mass1)
    merge_flag = np.zeros((numpairs),dtype="bool") # flag indicating whether the pair merges
    merge_snap = np.zeros((numpairs), dtype="int32") # the snapshot immediately after merger
    merge_redshift = np.zeros((numpairs), dtype="int32") # redshift immediately after merger
    infall_snap = np.zeros((numpairs), dtype="int32") # snapshot of first infall
    infall_redshift = np.zeros((numpairs)) # redshift of first infall
    pairkey = [] # a unique pair identifier 

    # 2D (per pair per snapshot):
    seps = np.zeros((numpairs, len(full_snaps))) # phys. separation at each snapshot
    seps_comov = np.zeros((numpairs, len(full_snaps))) # comoving separation
    seps_scaled = np.zeros((numpairs, len(full_snaps))) # scaled separation
    vels = np.zeros((numpairs, len(full_snaps))) # phys. velocity at each snapshot
    group_flag = np.zeros((numpairs,len(full_snaps)),dtype="bool") # True when the group numbers of prim & seco are the same
    rvir = np.zeros((numpairs, len(full_snaps))) # virial radius of the FoF group at each snapshot
    pos1, pos2 = np.zeros((2, numpairs, len(full_snaps), 3)) # x,y,z of each pair at each snapshot
    vel1, vel2 = np.zeros((2, numpairs, len(full_snaps), 3)) # vx,vy,vz of each pair at each snapshot 

    # fill seps and vels with np.nan, will be replaced with values after computation IF they can be computed
    seps.fill(np.NaN)
    seps_comov.fill(np.NaN)
    seps_scaled.fill(np.NaN)
    vels.fill(np.NaN)

    # read in the merger trees of the primary and secondary for each pair
    for ind in range(numpairs):
        # ----------------------------------------------------
        # get info and trees for primary and secondary
        primary_id = pairs['Sub1 ID'][ind]
        secondary_id = pairs['Sub2 ID'][ind]

        id1[ind] = primary_id
        id2[ind] = secondary_id
        mass1[ind] = pairs['Sub1 Mass'][ind]
        mass2[ind] = pairs['Sub2 Mass'][ind]
        stell1[ind] = pairs['Sub1 Stellar Mass'][ind]
        stell2[ind] = pairs['Sub2 Stellar Mass'][ind]
        smratio[ind] = pairs['Stellar Mass Ratio'][ind]
        groupnum[ind] = pairs['Group ID'][ind]

        treedict = {}
        for sub in [primary_id,secondary_id]:
            treedict[sub] = TraceMergerTree(snapshot=snapshot,physics="hydro",sim="TNG",subfindID=sub)

        tree_primary = treedict[primary_id].mergedbranch
        tree_secondary = treedict[secondary_id].mergedbranch

        # ----------------------------------------------------
        # check to see if halos merge
        root1 = tree_primary['RootDescendantID'][0]
        root2 = tree_secondary['RootDescendantID'][0]
        check_root = root1 == root2

        snap_mask1 = np.isin(tree_primary['SnapNum'],tree_secondary['SnapNum'])
        snap_mask2 = np.isin(tree_secondary['SnapNum'],tree_primary['SnapNum'])

        
        if check_root: # the primary and secondary DO merge in the simulation
            merge_flag[ind] = True

            # ----------------------------------------------------
            # identify the snapshot at which the halos have merged
            desc_mask1 = np.isin(tree_primary['DescendantID'],tree_secondary['DescendantID'])
            desc_mask2 = np.isin(tree_secondary['DescendantID'],tree_primary['DescendantID'])

            desc_overlap1 = tree_primary["SnapNum"][desc_mask1]
            desc_overlap2 = tree_secondary["SnapNum"][desc_mask2]

            if desc_overlap1[-1] == desc_overlap2[-1]:
                # this is the case when there are no snapshots skipped
                calc_snap = desc_overlap1[-1]+1
                merge_snap[ind] = calc_snap
                merge_redshift[ind] = convert['Redshift'][int(merge_snap[ind])]


            elif abs(desc_overlap1[-1]-desc_overlap2[-1])==1:
                # this is the case when the secondary snapshot temporarily enters 
                # the primary halo thus skipping a snapshot
                calc_snap = np.min(np.intersect1d(desc_overlap1,desc_overlap2))
                merge_snap[ind] = calc_snap
                merge_redshift[ind] = convert['Redshift'][int(merge_snap[ind])]

                
        else: # the primary and secondary do not merge in the simulation
            merge_snap[ind] = -1
            merge_redshift[ind] = -1


        # ----------------------------------------------------
        # collect separations, vels, etc. as function of snapshot
        for snap in full_snaps:
            scale = convert['Scale'][snap]

            # ----------------------------------------------------
            # test if both subhalos have data at the snapshot
            if (snap in tree_primary['SnapNum']) and (snap in tree_secondary['SnapNum']):
                # location of snapshot in trees
                loc1 = np.where(tree_primary['SnapNum'] == snap)[0]
                loc2 = np.where(tree_secondary['SnapNum'] == snap)[0]

                # ----------------------------------------------------
                # change flag to True if subhalos in same group
                grnum1 = tree_primary['SubhaloGrNr'][loc1]
                grnum2 = tree_secondary['SubhaloGrNr'][loc2]
                group_flag[ind][snap] = grnum1 == grnum2

                # ----------------------------------------------------
                # calculate separation between both halos
                subpos1 = tree_primary['SubhaloPos'][loc1] # comov
                subpos2 = tree_secondary['SubhaloPos'][loc2] # comov
                comoving_dist = np.linalg.norm(vector(subpos1,subpos2,75000))
                seps_comov[ind][snap] = comoving_dist
                seps[ind][snap]=comoving_dist*(scale)/little_h
                pos1[ind][snap] = subpos1
                pos2[ind][snap] = subpos2
                
                # ----------------------------------------------------
                # calculate the "scaled separation" compared to group r (note groups may be different)
                rvir_comov = tree_primary['Group_R_TopHat200'][loc1]
                rvir[ind][snap] = rvir_comov*scale/little_h
                seps_scaled[ind][snap] = comoving_dist/rvir_comov

                # ----------------------------------------------------
                # calculate velocity between both halos
                vel1[ind][snap] = tree_primary['SubhaloVel'][loc1]
                vel2[ind][snap] = tree_secondary['SubhaloVel'][loc2]
                rel_vel = np.linalg.norm(vel1[ind][snap]-vel2[ind][snap])
                vels[ind][snap]=rel_vel
                
        infall_snap[ind] = np.where(group_flag[ind]==True)[0][0]
        infall_redshift[ind] = convert['Redshift'][infall_snap[ind]]


        # ----------------------------------------------------
        # combine the last 6 digits of each of the first subhaloID from the primary and secondary tree
        pk_string = (str(tree_primary['SubhaloID'][-1])+str(tree_secondary['SubhaloID'][-1]))
        pairkey.append(pk_string.encode("utf-8"))
    
    # ----------------------------------------------------
    # make data structure and save to hdf5
    collection = {"Redshift":convert['Redshift'],
                  "Scale":convert['Scale'],
                  "Snapshot":convert['Snapshot'],
                  "GroupNum":groupnum,
                  "SubfindID1":id1,
                  "SubfindID2":id2,
                  "SubhaloMass1":mass1,
                  "SubhaloMass2":mass2,
                  "StellarMass1":stell1,
                  "StellarMass2":stell2,
                  "StellarMassRatio":smratio,
                  "MergeFlag":merge_flag, 
                  "MergeRedshift":merge_redshift,
                  "MergeSnapshot":merge_snap,
                  "InfallRedshift":infall_redshift,
                  "InfallSnapshot":infall_snap,
                  "PairKey":pairkey,
                  "GroupFlag":group_flag,
                  "GroupRvir":rvir,
                  "Separations":seps,
                  "SeparationsComoving":seps_comov,
                  "SeparationsScaled":seps_scaled,
                  "RelativeVelocity":vels,
                  "SubhaloPos1":pos1,
                  "SubhaloPos2":pos2,
                  "SubhaloVel1":vel1,
                  "SubhaloVel2":vel2}
    
    f = h5py.File(f"{paths.path_data}orbits/{masstype}mass_{pairtype}_{snapshot}.hdf5", 'w')

    info_dict = {"Redshift":"Redshift of snapshot",
                  "Scale":"Scale of snapshot",
                  "Snapshot":"Snapshot number",
                  "GroupNum":"FoF group number at snapshot",
                  "SubfindID1":"Subhalo ID of primary at selected redshift",
                  "SubfindID2":"Subhalo ID of secondary at selected redshift",
                  "SubhaloMass1":"Subhalo mass at selected redshift in 1e10*Msun",
                  "SubhaloMass2":"Subhalo mass at selected redshift in 1e10*Msun",
                  "StellarMass1":"Stellar mass from median AM at selected redshift",
                  "StellarMass2":"Stellar mass from median AM at selected redshift",
                  "StellarMassRatio":"Stellar mass ratio at selected redshift",
                  "MergeFlag":"True if subhalos will merge", 
                  "MergeRedshift":"Redshift of 'merger'",
                  "MergeSnapshot":"Snapshot at which 'merger' has occured",
                  "InfallRedshift":"First redshift where Group is the same",
                  "InfallSnapshot":"First snapshot where Group is the same",
                  "PairKey":"Unique identifying key for each pair (same between snapshots for same pair)",
                  "GroupFlag":"True if subhalos are in the same group",
                  "GroupRvir":"Radius of primary group in kpc",
                  "Separations":"Physical separation between pair in kpc",
                  "SeparationsComoving":"Comoving separation between pair in ckpc/h",
                  "SeparationsScaled":"Separation scaled by radius of primary group, dimensionless",
                  "RelativeVelocity":"Relative velocity between pair in km/s",
                  "SubhaloPos1":"Position of subhalo in ckpc/h",
                  "SubhaloPos2":"Position of subhalo in ckpc/h",
                  "SubhaloVel1":"Velocity of subhalo in km/s",
                  "SubhaloVel2":"Velocity of subhalo in km/s"
                }

    for key, val in collection.items():

        if key=='PairKey':
            dset = f.create_dataset(f'/{key}',
                                    shape=(len(pairkey),),
                                    dtype=h5py.string_dtype(encoding='utf-8', length=None))
            dset.attrs[key] = info_dict[key]
            dset[:]=np.array(pairkey)
        else:
            print(key)
            valv = np.array(val)
            dset = f.create_dataset(f'/{key}', 
                                shape=valv.shape,
                                dtype=valv.dtype)
            dset.attrs[key] = info_dict[key]
            dset[:] = valv

    f.close()    
    print("Saved orbit data for snapshot ",snapshot)
    