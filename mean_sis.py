'''
Created on Feb 18, 2014

@author: BEC3
'''
from readsis import loadimg, write_raw_image
import numpy, settings, numpy.ma
from __builtin__ import len


import os

def get_filepaths(directory):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths  # Self-explanatory.


def create_reference(sis_to_mean):
    meanpath=settings.referencefile
    
    for f in sis_to_mean[:]:
        if not f.endswith(".sis"):
            sis_to_mean.remove(f)
            
    if len(sis_to_mean)>0:
        tot_K=None
        tot_Na=None
        
        for file in sis_to_mean:
            K,Na=loadimg(file)
            if tot_K!=None:
                tot_K+=K
            else:
                tot_K=K
            if tot_Na!=None:
                tot_Na+=Na
            else:
                tot_Na=Na
        
        tot_K/=len (sis_to_mean)
        tot_Na/=len( sis_to_mean)
        
        tot=numpy.concatenate((numpy.ma.filled(tot_K,0), numpy.ma.filled(tot_Na,0)))
        tot=(tot+1)*6553.6
        
        write_raw_image(meanpath,tot)

if __name__ == '__main__':
    
    
    path="D:/SIScam/SIScamProgram/Prog/img/2014/2014-02-18/bitmaps"
    
    sis_to_mean=get_filepaths(path)
    create_reference(sis_to_mean)