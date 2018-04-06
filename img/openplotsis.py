#!/usr/bin/python3
# -*- coding: utf-8 -*-# -*- coding: utf-8 -*-
import numpy as np

class RawSis():
    
    def __init__(self, filename):
        
        im0, im1, im, raw = self.readsis(filename)
        
        self.im0 = im0
        self.im1 = im1
        self.im_full  = im
        self.raw = raw
#        self.height, self.width = hw
        
    
    def readsis(self,filename):
        f = open(filename, 'rb')  #apre in binario
        rawdata = np.fromfile(f,'H').astype(int)
        f.close
        
        width=rawdata[6]  # N cols
        height=rawdata[5] # N rows
        #rispetto ad octave, gli indici cambiano (python Ã¨ 0-based)
    
        image = rawdata[-width*height : ]
        image.resize(height,width)
               
        #in matlab, poichÃ© reshape legge per colonne, era necessario scambiare gli indici e poi trasporre la matrice. Ora non piÃ¹
        
        im0 = image[:height//2, :]
        im1 = image[height//2:, :]
        
        
        return im0, im1, image, rawdata #, image.shape
 
    
if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    sis = RawSis(sys.argv[1])    
    fig, (ax1, ax2) = plt.subplots(2,1)
    dic = {'cmap': 'gist_stern'}
    ax1.imshow(sis.im0, **dic)
    ax2.imshow(sis.im1, **dic)

    fig.show()
       
        
    
    
    
    