#!/usr/bin/env python
"""Routines to determine and manipulate the file name formats of various
output files."""

__author__ = "Mark Vogelsberger, Paul Torrey and contributing authors"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Mark Vogelsberger, Paul Torrey and contributing authors"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@mit.harvard.edu"
__status__ = "Beta -- forever."

import glob
import os
import sys
import numpy as np

def get_snap_filenames(base, snap_prefix, num):
    """Return a list of paths (without file extensions) to snapshot files
    corresponding to the snapshot indicated by num.  The list may consist of
    one or multiple files, depending on how the snapshot was written.
    """
    num_pad = str(num).zfill(3)
    path = "%s/%s_%s" % (base, snap_prefix, num_pad)
    if os.path.exists("%s.hdf5" % path):
        return [path]

    path = "%s/snapdir_%s/%s_%s" % (base, num_pad, snap_prefix, num_pad)
    if os.path.exists("%s.0.hdf5" % path):
        files = np.array(glob.glob("%s.*.hdf5" % path))
        basenames = [os.path.basename(this_file) for this_file in files]
        filenrs =  np.array([int(this_file[this_file.index('.')+1:this_file.index('.hdf5')]) for this_file in basenames])
        files = files[ np.argsort(filenrs) ]
        return [os.path.splitext(x)[0] for x in files]

    raise RuntimeError("Cannot find snapshot number %s in %s." % (num, base))



#def return_subfind_filebase(basedir, snapnum, name, filenum):
#    filebase = basedir + "/groups_" + str(snapnum).zfill(3) + "/" + name + "_" + str(snapnum).zfill(3) + "."
#    curfile = filebase + str(filenum) + ".hdf5"
#    
#    if not os.path.exists(curfile):
#        filebase = basedir + "/" + name + "_" + str(snapnum).zfill(3)
#        curfile = filebase + ".hdf5"
#    if not os.path.exists(curfile):
#        filebase = basedir + "/output/" + name + "_" + str(snapnum).zfill(3)
#        curfile = filebase + ".hdf5"
#    if not os.path.exists(curfile):
#        filebase = basedir + "output/groups_" + str(snapnum).zfill(3) + "/" + name + "_" + str(snapnum).zfill(3) + "."
#        curfile = filebase + str(filenum) + ".hdf5"
#    
#    if not os.path.exists(curfile):
#        print "file not found:", curfile
#        sys.exit()
#    
#    return filebase, curfile

def return_subfind_filebase(basedir, snapnum, name, filenum):
    return return_general_filebase( basedir, snapnum, name, filenum, '/groups_')

def return_snap_filebase(basedir, snapnum, name="snap", filenum=0):
    return return_general_filebase( basedir, snapnum, name, filenum, '/snapdir_')

def return_general_filebase(basedir, snapnum, name, filenum, string):
    filebase = basedir + string + str(snapnum).zfill(3) + "/" + name + "_" + str(snapnum).zfill(3)
    curfile = filebase + "." + str(filenum) + ".hdf5"

    if not os.path.exists(curfile):
        filebase = basedir + "/" + name + "_" + str(snapnum).zfill(3)
        curfile = filebase + ".hdf5"
    if not os.path.exists(curfile):
        filebase = basedir + "/output/" + name + "_" + str(snapnum).zfill(3)
        curfile = filebase + ".hdf5"
    if not os.path.exists(curfile):
        filebase = basedir + "output" + string + str(snapnum).zfill(3) + "/" + name + "_" + str(snapnum).zfill(3) 
        curfile = filebase + "." + str(filenum) + ".hdf5"
    if not os.path.exists(curfile):
        print("file not found:", curfile)
        sys.exit()

    return filebase, curfile


def return_galprop_file(basedir, snapnum):
    snaptag='000'+str(snapnum)
    snaptag=snaptag[-3:]
    file=basedir+'/postprocessing/galprop/galprop_'+snaptag+'.hdf5'
    if os.path.exists(file):
        return file
    file=basedir+'/galprop/galprop_'+snaptag+'.hdf5'
    if os.path.exists(file):
        return file
    file=basedir+'/galprop_'+snaptag+'.hdf5'
    if os.path.exists(file):
        return file
    print("Galprop File Not Found... (see util.naming for attempted paths)")
    return None


def partTypeNum(partType):
    """ Mapping between common names and numeric particle types. """
    if str(partType).isdigit():
        return int(partType)
    
    if str(partType).lower() in ['gas','cells']:
        return 0
    if str(partType).lower() in ['dm','darkmatter']:
        return 1
    if str(partType).lower() in ['tracer','tracers','tracermc','trmc']:
        return 3
    if str(partType).lower() in ['star','stars','stellar']:
        return 4 # only those with GFM_StellarFormationTime>0
    if str(partType).lower() in ['wind']:
        return 4 # only those with GFM_StellarFormationTime<0
    if str(partType).lower() in ['bh','bhs','blackhole','blackholes']:
        return 5
    
    raise Exception("Unknown particle type name.")


