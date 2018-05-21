import os

#basedir = r"c:/fit/siscam"
#basedir = os.getcwd()
basedir = r"c:/SIScam/SIScamProgram/Prog"

#full path to image file
imgfold = r"c:/SIScam/SIScamProgram/Prog/img"
imagefile = os.path.join(imgfold, "test.sis")
watchedfiles = [os.path.join(imgfold, f) for f in ['test_0.sis', 'test_1.sis']]

rawimage1file = imagefile
rawimage2file = imagefile
rawimage3file = imagefile

referencefile = r"c:/SIScam/SIScamProgram/Prog/img/reference.sis"

#where to save images

imagesavepath = os.path.normpath("d:/SIScam/SIScamProgram/Prog/img")

#icons, etc.
bitmappath = os.path.join(basedir, 'bitmaps')

#directory to store template files
templatedir = os.path.join(basedir, 'templates')

##acquire
useTheta = True
useBluefox = True
useSony = False

#settings for Theta Systems cam

configfile = r'c:/SIScam/SIScamProgram/Prog/config.ini'

#usePseudoCam = False
usePseudoCam = True
