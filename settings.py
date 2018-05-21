import os

#basedir = r"c:/fit/siscam"
#basedir = os.getcwd()
basedir = r"c:/SIScam/SIScamProgram/Prog"

#full path to image file
#imagefile = r"c:/fit/ocf/fit.sis"
#imagefile = r"./img/test.sis"
imagefile = r"c:/SIScam/SIScamProgram/Prog/img/test.sis"
#imagefile = r"./img/fit.sis"
#imagefile = r"./img/20100319-LBside-0010.sis"
#rawimage1file = r"c:/fit/ocf/rawimg1.sis"
#rawimage2file = r"c:/fit/ocf/rawimg2.sis"
#rawimage3file = r"c:/fit/ocf/rawimg3.sis"
rawimage1file = imagefile
rawimage2file = imagefile
rawimage3file = imagefile

referencefile = r"c:/SIScam/SIScamProgram/Prog/img/reference.sis"

#where to save images
#imagesavepath = r"./img"
#imagesavepath = r"d:/SIScam/SIScamProgram/Prog/img"
imagesavepath = os.path.normpath("d:/SIScam/SIScamProgram/Prog/img")

#icons, etc.
bitmappath = os.path.join(basedir, 'bitmaps')

#directory to store template files
templatedir = os.path.join(basedir, 'templates')

##acquire
useTheta = True
#useTheta = False
useBluefox = True
#useBluefox = False
useSony = False

#settings for Theta Systems cam
#configfile = r"c:\WinSIS6\py\config.ini"
#configfile = r"./WinSIS6/py/config.ini"
configfile = r'c:/SIScam/SIScamProgram/Prog/config.ini'

#usePseudoCam = False
usePseudoCam = True
