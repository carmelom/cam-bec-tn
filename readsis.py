#!/usr/bin/python
#-*- coding: latin-1 -*-
"""Load image files which have been saved by WinSIS."""

from __future__ import with_statement
import numpy
import numpy as np
from numpy import fromfile, uint16, log, fromstring, ma
import os.path, sys, time
import settings_readsis

def read(filename):
    '''
    INTRODUCED ON 2016-10-19 to read from ELENA+Hamamatsu
    
    Read sis files, both old version and SisV2
    Parameters
    ----------
    filename : stringa con il nome o il path relativo del file.
    Returns im1, im2, image, rawdata, stringa, block
    -------
    im1 : ndarray 2D
        first half of the image [righe 0 : height/2-1].
    im2 : ndarray 2D
        second half of the image [righe height/2 : height-1].
    image : ndarray 2D
        whole image
    rawdata : ndarray 1D
        raw data read from file
    stringa : string
        comment and datestamp of the image
    block : tuple (Bheight,Bwidth)
        Bheight : int
            y dimension of the sub-block
        Bwidth : int
            x dimension of the sub-block
    Notes
    -----
    NB all the outputs are slices of the raw data: any modifications of them will be reflected also in the linked rawdata elements
    '''
    f = open(filename, 'rb')        # open in reading and binary mode
    header = f.read(512)            # read the fixed size header: 512 bytes
    #if header[0] == 48:             # 48 is ASCII code for '0': all older sis start with 0
    sisVersion = '  sis1'
    datestamp = ' we do not know'
    commitProg = '  nothing'
    comment = ' nothing'
    # elif header[0:5] == b'SisV2':   # SisV2 print header
        # sisVersion = str(header[0:5])
        # datestamp = str(header[18:37])
        # commitProg = str(header[40:48])
        # comment = str(header[48:])
    # else:
        # print('What fuck are you opening?')
        # return True

    print("You are opening the " + sisVersion[2:-1] + " file: " + filename)
    f.close()
    # Sometimes it gives error if the file is opened a sigle time

    f = open(filename, 'rb')                        # open in reading and binary mode
    rawdata = np.fromfile(f,'H').astype(int)
    # 'H' = uint16
    # types are listed in np.typeDict
    # put in an array the data formatted uint16 and casted to int
    ''' NB fundamental to cast to int:
        unsigned short gives overflow '''
    f.close()

    # Dimension of the whole image
    width=rawdata[6]  # N cols
    height=rawdata[5] # N rows

    # Dimension of the sub-blocks
    Bwidth=rawdata[8]  # N subcols
    Bheight=rawdata[7] # N subrows

    # Length of comments
    ls = rawdata[19]

    # Reading the images
    # image = rawdata[-width*height:]
    # image.resize(height,width)

    image = np.resize(rawdata[-width*height:], (height,width))
    return image.astype(np.float32)

def OLDread(filename):
    fid = file(filename, 'rb')
    foo = fid.read(10)
    height = int(fromfile(fid, dtype = uint16, count = 1))
    width  = int(fromfile(fid, dtype = uint16, count = 1))
    fid.read(182)
    xoff = int(fromfile(fid, dtype = uint16, count = 1))
    yoff = int(fromfile(fid, dtype = uint16, count = 1))

    #img = fromfile(fid, dtype = uint16, count = width*height)
    data = read_fid_full(fid, size = width*height*2, timeout = 5)
    img = fromstring(data, dtype = uint16, count = width*height)
    
    img.shape = (height, width)

    fid.close()
    #return img.astype(numpy.float_)
    return img.astype(numpy.float32)

def read_fid_full(fid, size, timeout = 1):
    numbytesread = 0
    result = ''
    starttime = time.clock()
    while numbytesread < size:
        result += fid.read(size - numbytesread)
        numread = len(result)

        if time.clock() - starttime > timeout or numread>=size:
            break
        time.sleep(0.1)
    return result

def write_raw_image(filename, img):
    fid = file(filename, 'wb')
    #fid.write(' '*10)   #############original
    fid.write('0'*10)   #############01-02-2013
    height, width = img.shape
    ha = numpy.array([height], dtype = numpy.uint16)
    ha.tofile(fid)
    wa = numpy.array([width], dtype = numpy.uint16)
    wa.tofile(fid)
    #fid.write(' '*182)   #############original
    fid.write('0'*(182))   #############01-02-2013
    #fid.write(' '*4) #TODO: xoff, yoff   #############original
    fid.write('0'*4)   #############01-02-2013
    if img.dtype == numpy.uint16:
        img.tofile(fid)
    else:
        img.astype(numpy.uint16).tofile(fid)
    fid.close()

def loadimg(filename):
    img = read(filename)
    
    h, w = img.shape
    
    imgNa = img[:h/2]
    imgK = img[h/2:]
    #return imgNa/6553.6 - 1.0, imgK/6553.6 - 1.0
    return ma.masked_where(imgNa==0, imgNa/6553.6-1.0), \
           ma.masked_where(imgK==0, imgK/6553.6-1.0)

def loadimg3(path):
    img1 = read(os.path.join(path, 'PIC1.SIS'))
    img2 = read(os.path.join(path, 'PIC2.SIS'))
    img3 = read(os.path.join(path, 'PIC3.SIS'))

    img = - (log(img1 - img3) - log(img2 - img3))
    return img[:1040], img[1040:], 

def test_write_read():
    #imgNa, imgNa = loadimg('img/test.sis')   #####added 14-12-2012
    #imgK, imgK = loadimg('img/test_1.sis')   #####added 14-12-2012
    imgK, imgNa = loadimg('img/test.sis')   #####original 14-12-2012
    #imgK, imgNa = loadimg3('this string is for nothing')

    rawimg = (6553.6*(imgK + 1)).astype(numpy.uint16)   # 1000
    
    write_raw_image('img/testsave.sis', rawimg)
    rawimgsaved = read('img/testsave.sis')

def test_read(filename):
    print "loading..."
    sys.stdout.flush()
    loadimg(filename)
    print "done"
    sys.stdout.flush()

def test_save(filename, img):
    print "saving..."
    sys.stdout.flush()
    write_raw_image(filename, img)
    print "done"
    sys.stdout.flush()

def simultanous_write_read():
    import threading
    imgK, imgNa = loadimg(settings_readsis.imagefile)
    rawimg = (6553.6*(imgK + 1)).astype(numpy.uint16)  # 1000

    savethread = threading.Thread(target = test_save,
                                  args = ('./img/testsave.sis', rawimg),
                                  )

    loadthread = threading.Thread(target = test_read,
                                  args = ('./img/testsave.sis',),
                                  )

    savethread.run()
    loadthread.run()


if __name__ == '__main__':
    #test_write_read()
    simultanous_write_read()
    
