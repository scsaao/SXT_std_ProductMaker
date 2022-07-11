
    

def lcplotter(lcf, binsize = 60, en = "0p3to7p0", outfile = "lcout.pdf", LCOutTable = 'lc.ipac',
				figsize = [10,9]) :
	import numpy as np
	import os
	from astropy.io import fits
	from astropy.table import Table, Column
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter
	from astropy.io import ascii
	from astropy.io.ascii import masked
	
	f = fits.open(lcf, ignore_missing_end=True)
	lcdt = f[1].data
	lchdr = f[0].header
	lchdr2 = f[1].header
	f.close()

	tstart = lchdr['TSTART']
	tstop = lchdr['TSTOP']
	mjdref = lchdr['MJDREFI']
	srcname = lchdr['object']
	exposure = lchdr['exposure']
	ra, dec = lchdr2['RA_NOM'], lchdr2['DEC_NOM']
	mjdstart = tstart/86400. + mjdref
	obsid = lchdr2['OBS_ID']
	binsize = lchdr2['TIMEDEL']
	pnt_coord = (ra, dec)
	E_Range = en
	target = srcname						
	commentstr = f"# NOMINAL POINTING : {(ra, dec)}\n"
	commentstr += f"#MJDSTART : {mjdstart}\n"
	commentstr += f"#EnergyString : {en}\n"

	#fill = [(masked, 'N/A', 'ra'), (masked, -999, 'sptype')]

	print(commentstr)
	enflag = en.split("to")[0].split("p")[0] +"."+ en.split("to")[0].split("p")[1] + "-" + en.split("to")[1].split("p")[0] + "." +en.split("to")[1].split("p")[1]
	#getting source name without white space
	if len(srcname.split(" ")) > 1 :
		sname = srcname.split(" ")[0]
		for u in range(len(srcname.split(" "))-1) :
			sname = "{}{}".format(sname,srcname.split(" ")[u+1])
		srcname = sname
		
	
	if len(lcdt) > 3 :
		lcdt = lcdt[lcdt['FRACEXP'] >= 0.9]
		lcdtflxm = np.nanmean(lcdt['RATE']); lcdtflsd = np.nanstd(lcdt['RATE'])
		lcdt = lcdt[lcdt['RATE'] >= (lcdtflxm - 3.5*lcdtflsd)]
		lcdt = lcdt[lcdt['RATE'] <= (lcdtflxm + 3.5*lcdtflsd)]

		enstring = en.split('to')[0].split('p')[0] + '.' + en.split('to')[0].split('p')[1] +'-'+ en.split('to')[1].split('p')[0] + '.' + en.split('to')[1].split('p')[1]

		#...writing to file----------------- 
		tblout = Table()
		tblout.add_column(Column((lcdt['TIME']/86400.) + mjdstart , name = "TIME"))
		tblout.add_column(Column(lcdt['RATE'], name = "RATE"))
		tblout.add_column(Column(lcdt['ERROR'], name = "ERROR"))
		tblout.add_column(Column(lcdt['FRACEXP'], name = "FRACEXP"))
		tblout['TIME'].unit = 'MJD' 
		tblout['RATE'].unit = 'c/s' 
		tblout['ERROR'].unit = 'c/s'
		
		keydict = {}
		for keywrds in ['target', 'E_Range', 'obsid', 'tstart', 'tstop', 'exposure', 
						'mjdstart', 'binsize', 'pnt_coord']:
			keydict[keywrds.upper()] = {'value': eval(keywrds)}
		DictMeta = {"comments": "SXT Light Curve", 'keywords' : keydict}	
		tblout.meta = DictMeta
		 
		tblout.write(LCOutTable, format="ascii.ipac", comment = commentstr)
		#tblout.close()
		#-----------------------------------

		fig, ax = plt.subplots(1, sharex=True, figsize = figsize)
		x = lcdt['TIME']; y = lcdt['RATE']; ey = lcdt['ERROR']

		xlab = "Time [s; MET] - {}".format(round(mjdstart,1))
		ylab = "Rate [c/s]"

		ax.errorbar(x, y, xerr= binsize/2, yerr= ey, fmt="*k", markersize = 8, 
						label = "binsize : {} s".format(binsize))

		ylim = [y.min() - 0.3*(y.max() - y.min()),  y.max() + 0.3*(y.max() - y.min())]
		xlim = [x.min() - 0.1*(x.max() - x.min()),  x.max() + 0.2*(x.max() - x.min())]
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		ax.set_xlabel(xlab, fontsize = 23)
		ax.set_ylabel(ylab, fontsize = 23)
		text2use4t = f"Object:{srcname}; E:{enflag} keV; Exp.:{round(exposure/1e3,2)} ks"
		ax.set_title(f"SXT Lightcurve; {text2use4t}", fontsize = 21)
		major_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/7.))
		minor_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/35.))
		ax.tick_params(axis = 'both', which = 'major', labelsize = 19)
		ax.tick_params(axis = 'both', which = 'minor', labelsize = 17)
		ax.set_xticks(np.round(major_ticks,1))
		ax.set_xticks(minor_ticks, minor = True)
		ax.grid(which = 'major', alpha = 0.4, linewidth=1, color = 'k', linestyle = '--')
		ax.grid(which = 'minor', alpha = 0.3, linewidth=1, color = 'k', linestyle = '-.')
		ax.legend(loc=1)
		plt.savefig(outfile)
		fig.clf()
		plt.close()

		if os.path.exists(outfile) :
			status = "LCPlotter succeded in making lighcurve plots with name : \n {}".format(outfile)
		else :
			status = "Warning:: Something fissy happend could not write the curve plot"

	else :
		status = "Error:: The length of the lc data is significantly lower (=< 3), so skipping the plotting part"
	#plt.show()
	return status, ax


def phaplotter(phafname, grpmin = 60, outfile = "lcout.pdf", sxtrspdir = "./", 
				respfilename = "sxt_pc_mat_g0to12.rmf", ancrfilename = "sxt_onaxis_rad05.arf", 
				backfilename = "All6phaadded_spec.pha", sname = "None", figsize=[10,9],
				phatyp = 'unbinned', bkgflag = False) :
	import numpy as np
	import os
	from astropy.io import fits
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter

	backfilename = os.path.join(sxtrspdir, backfilename)
	respfilename = os.path.join(sxtrspdir, respfilename)
	ancrfilename = os.path.join(sxtrspdir, ancrfilename)
    
	if phatyp in ['g', 'grpd', 'grouped', 'grppha']:
		grpflag = True
	else:
		grpflag = False
	
	if grpflag:
		grpdphaname = phafname.split('.pha')[0]
		grpdphaname = f"{grpdphaname}_gr.pha"
		cmd = f"grppha {phafname} "
		cmd += f" COMM='group min {grpmin} &"
		if bkgflag: 
			cmd += f" chkey backfile {backfilename} &"
		cmd += f" chkey respfile {respfilename} &"
		cmd += f"chkey ancrfile {ancrfilename} & exit"
		os.system(cmd)
	else:
		grpdphaname = phafname 
	#os.system("grppha {} {}_gr.pha COMM='group min {} & chkey respfile {} & chkey ancrfile {}"
	#" & exit'".format(phafname, phafname.split('.pha')[0]"
	#", grpmin, os.path.join(sxtrspdir, respfilename), os.path.join(sxtrspdir, ancrfilename)))"

	grpphaf = fits.open(grpdphaname, ignore_missing_end=True)
	grpphafhdr = grpphaf[0].header
    
	try:
		objname = grpphafhdr['OBJECT']
	except:
		objname = sname
		
	grpphafdata = grpphaf[1].data
	grpphaf.close()
	#f = fits.open('phaf')
	tstart = grpphafhdr['TSTART']
	tstop = grpphafhdr['TSTOP']
	mjdref = grpphafhdr['MJDREFI']
	mjdstart = tstart/86400. + mjdref
	dateobs = grpphafhdr['date-obs']
	timeobs = grpphafhdr['time-obs']
	timeobs = "{}:{}:{}".format(timeobs.split(":")[0], timeobs.split(":")[1], 
									int(float(timeobs.split(":")[2])))
	exp = grpphafhdr['EXPOSURE']
	obsid = grpphafhdr['OBS_ID']
	exp2prnt = str(round(exp/1e3,1)).split(".")[0]+"p"+str(round(exp/1e3,1)).split(".")[1]+"ks"
    
	#making tcl script if not in dir...
    
	if os.path.exists("tcscrpt2run.tcl") :
		os.remove("tcscrpt2run.tcl")

	tclf = open("tcscrpt2run.tcl",'w')
	tclf.write("cpd /xw \n")
	tclf.write(f"data {grpdphaname} \n")
	if grpflag:
		
		tclf.write("setpl e \n")
		tclf.write("ig **-0.3 7.0-** \n")
		tclf.write("ig bad \n")
		XLABEL = "Energy [keV]"	
	else:
		XLABEL = "Energy [Channels]"	
		
	tclf.write("plot ldata \n")
	tclf.write("setplot delete all \n")
	tclf.write("setplot device /xw \n")
	tclf.write(f"setplot command we tcscrpt2runOut \n")
	tclf.write("plot ldata")
	tclf.close()

	# Running tcl scripts for generating qdp for spectrum
	os.system('xspec tcscrpt2run.tcl')
	os.remove('tcscrpt2runOut.pco')
	os.remove('tcscrpt2run.tcl')
	phadt = np.loadtxt('tcscrpt2runOut.qdp', skiprows=3)
	os.remove('tcscrpt2runOut.qdp')
	if len(phadt) > 3 :
		fig, ax = plt.subplots(1, sharex=True, figsize = figsize)
		x = phadt[:,0] ; ex = phadt[:,1]; y = phadt[:,2]; ey = phadt[:,3]
    
		xlab = XLABEL
		ylab = r"Norm. Counts s$^{-1}$ keV$^{-1}$"

		ax.errorbar(x, y, xerr= ex, yerr= ey, fmt="*k", 
						label = f"Date-OBS : {dateobs}T{timeobs} \n Exp. : {np.round(exp/1e3,1)} ks",
						markersize = 8)
		#print (np.column_stack((x,y)))
		if (y.min()/y.max()) > 1e-3 :
			ylim0 = y.min()
		else :
			ylim0 = 1e-4

		ylim = [ylim0,  y.max() + 0.4*(y.max() - y.min())]
		if grpflag:
			xlim = [0.25,  8.0]
		else:
			xlim = [25, 800]
			ylim = [ylim0,  2e-2]	
			ylim = [ylim0,  3e-3]
		ax.set_xlim(xlim)
		ax.set_ylim(ylim)
		ax.set_xlabel(xlab, fontsize = 17)
		ax.set_ylabel(ylab, fontsize = 17)
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.text(0.2,0.97, f"OBSID: {obsid}", horizontalalignment='center', 
					verticalalignment='center', transform=ax.transAxes)
					
		#major_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/5.))
		#minor_ticks = np.arange(xlim[0],xlim[1],int((xlim[1]-xlim[0])/25.))
		if bkgflag:
			ax.set_title(f"SXT Spectrum; [ Target: {objname} ]", fontsize = 19)
		else:
			ax.set_title(f"SXT Spectrum; No Bkg. Sub. [Target: {objname}]", fontsize = 19)	
		ax.tick_params(axis = 'both', which = 'major', labelsize = 17)
		ax.tick_params(axis = 'both', which = 'minor', labelsize = 15)
		#ax.set_xticks(np.round(major_ticks,1))
		#ax.set_xticks(minor_ticks, minor = True)
		ax.grid(which = 'major', alpha = 0.3, linewidth=1, color = 'k', linestyle = ':')
		ax.grid(which = 'minor', alpha = 0.1, linewidth=1, color = 'k', linestyle = ':')
		ax.legend(loc=3)
		plt.savefig(outfile)
		fig.clf()
		plt.close()

		if os.path.exists(outfile) :
			status = "PHAPlotter succeded in making spectral plot with name : \n {}".format(outfile)
		else :
			status = "Warning :: Something fissy happend could not write the spectral plot"
	else :
		status = "Error:: The length of the spectrum files is less than 3, so skipping the plot part..." 
	return status, ax
    #plt.show()


def StdProductMaker(evtfile,regdir="./",regfile=None,CHANMIN=30, CHANMAX=701, grade_flag="0-12",
						productstem="SrcNameOut", curvebinsize = 120, xyname = "RAW", 
						timefile = None, intensity = None):
						
						 
    
    import os
    
    #auxParDict 
    if regfile == None :
        print ("The product is contaminated with corner sources")
        status = "The product is contaminated with corner sources"
        regfile = ''
    else :
        status = "The product is made with your input"

    xco_outstr="{}_scrpt.xco".format(productstem)
    
    regfInp=os.path.join(regdir,regfile)
    
    evtdir="./"
    
    #Basic settings in .xco file
    #Calling XSELECT
    fl_xco=open(str(xco_outstr),'a')  
    fl_xco.write("xsel\n")
    
    #READING EVENTS FILE
    fl_xco.write("read events\n")
    fl_xco.write(f"{evtdir}\n")
    fl_xco.write(f"{evtfile}\n")
    fl_xco.write("yes\n")
    
    #SETTING UP THE XYNAME (RAW v/s SKY)
    if str(xyname).upper() == "RAW" :
        fl_xco.write("set xyname RAWX RAWY\n")
    else :
        fl_xco.write("set xyname X Y\n")
        
    #SETTING UP THE PHANAME (PI v/s PHA)
    fl_xco.write("set PHANAME PI\n")
    
    #CALLING EXTRACT FOR CHECK PRIMARY RUN
    fl_xco.write("extract all\n")
    
    #APPLYING CHANNEL FILTER
    fl_xco.write(f"filter pha_cutoff {CHANMIN} {CHANMAX} \n")
    
    #APPLY GRADE FILTER
    if grade_flag != None:
        fl_xco.write(f"filter grade {grade_flag}\n")
    
    #APPLY INTENSITY FILTER
    if intensity != None:
        fl_xco.write(f"filter intensity {intensity}\n")
    
    #APPLY TIME FILTER USING FILE (T1 T2)    
    if timefile != None:
	    fl_xco.write(f"filter time file {timefile}\n")
	
	#APPLY REGION FILTER	   
    fl_xco.write(f"filter region {regfInp} \n")
    
    #SAVE IMAGE    
    fl_xco.write("extract image\n")
    fl_xco.write(f"save image {productstem}.img \n")

    #SAVE LIGHTCURVE WITH BINSIZE = delta(T)
    fl_xco.write(f"set binsize {curvebinsize} \n")
    fl_xco.write("extract curve\n")
    fl_xco.write(f"save curve {productstem}.lc \n")
    fl_xco.write("clear pha_cutoff\n")
    
    #SAVE SPECTRUM
    fl_xco.write("extract spectrum\n")
    fl_xco.write(f"save spectrum {productstem}.pha \n")
    
    #QUITTING WITHOUT SAVING THE SESSION
    fl_xco.write("quit\n")
    fl_xco.write("no\n")
    fl_xco.close()
    
    os.system(('xselect @%s')%(xco_outstr))
    os.system("rm -rf {}".format(regfInp))
    return status 


def gkern(kernlen=21, nsig=3):
    import scipy.stats as st
    import numpy as np

    """Returns a 2D Gaussian kernel array."""
    interval = (2*nsig+1.)/(kernlen)
    x = np.linspace(-nsig-interval/2., nsig+interval/2., kernlen+1)
    kern1d = np.diff(st.norm.cdf(x))
    kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
    kernel = kernel_raw/kernel_raw.sum()
    return kernel

def imgaeplotter(infile = "imgname", outfile = "outfile.png", figsize = [10,9], krnrad = 15,
					krnsig = 3, imgflag = 's') :
    import matplotlib.pyplot as plt
    import os, seaborn
    from astropy.io import fits
    import scipy
    import scipy.signal
    import numpy as np
    import matplotlib.colors as colors
    
    imgf = fits.open(infile, ignore_missing_end=True)
    imgfdata = imgf[0].data
    imgfhdr = imgf[0].header
    imgf.close()
    
    mrgd_exp = imgfhdr['EXPOSURE']
    srcname = imgfhdr['object']
    obsid = imgfhdr['OBS_ID']
    
    if imgflag == 'r':
        ImgFrmSize = [0, 600]
        ylab = "Dec [DET; pix]"
        xlab = "RA [DET; pix]"
    else:
        ImgFrmSize = [0, 1000]
        ylab = "Dec [SKY; pix]"
        xlab = "RA [SKY; pix]"
			
    if len(srcname.split(" ")) > 1 :
        sname = srcname.split(" ")[0]
        for u in range(len(srcname.split(" "))-1) :
            sname = "{}{}".format(sname,srcname.split(" ")[u+1])
        srcname = sname

    r = scipy.signal.convolve2d(imgfdata, gkern(kernlen = krnrad, 
									nsig = krnsig), mode='same')
    fig, ax = plt.subplots(figsize = figsize)
    
    ax.set_title(f"OBSID: {obsid}, Exp: {str(round(mrgd_exp/1e3,1))} ks", fontsize = 20)
    img1 = ax.imshow(r, origin = 'lower', interpolation='nearest', cmap=plt.cm.jet,
						norm = colors.LogNorm())
				
    fig.colorbar(img1, ax=ax)
        			
    #sub.imshow(r, interpolation='nearest', vmin=.0001, vmax=r.max()-0.9*r.max() ,cmap=plt.cm.PuBuGn)
    
    ax.set_xlim(ImgFrmSize[0], ImgFrmSize[1])
    ax.set_ylim(ImgFrmSize[0], ImgFrmSize[1])
    ax.set_ylabel(ylab, fontsize = 23)
    ax.set_xlabel(xlab, fontsize = 23)
    
    #seaborn.set_style("ticks")
    #ax.text(0.1,0.9, f"OBSID: {obsid}", horizontalalignment='center', fontsize = 7,
	#			verticalalignment='center', transform=ax.transAxes)
    print ("The Device has made your plot...save it as per your choice...")
    plt.tick_params(axis='both',which='major',labelsize=21)
    plt.tick_params(axis='both',which='minor',labelsize=19)
    fig.savefig("{}".format(outfile))
    fig.clf()
    plt.close()

    if os.path.exists(outfile) :
        status = "The image is succesfully written as file: {}".format(outfile)
    else :
        status = "The image maker is failed "
    return status, ax 


def time_transform(mjd) :
    from astropy.time import Time
    mjdtime = Time(mjd, format='mjd')
    uttime = mjdtime.isot
    decimalyear = mjdtime.decimalyear
    decimalyear = [round(float(dy),2) for dy in decimalyear if dy]
    return  uttime,  decimalyear

def sxt_met2mjd(timevec) :
    mjdt = []
    for t in timevec :
        mjdt.append((t/86400.) + 55197.0)
    return mjdt
    
    
def GetCentralcoordLookup(evtfile, rawx=300.0, rawy = 300.0, coordtol=2e0, stepcoordtol=5e-1) :

	import numpy as np
	from astropy.io import fits
	from astropy.table import Table

	evtfl = fits.open(evtfile)
	evthdu = evtfl[1].data
	evtfl.close()

	evttbl = Table(evthdu)
	N = 0

	while N <= 3 :
		xrawlowlim, xrawuplim = rawx - coordtol, rawx + coordtol 
		yrawlowlim, yrawuplim = rawy - coordtol, rawy + coordtol
		fltrdhdu = evttbl[(evttbl['RAWX']>=xrawlowlim) & (evttbl['RAWX']<=xrawuplim) & 
							(evttbl['RAWY']>=yrawlowlim) & (evttbl['RAWY']<=yrawuplim)]

		coordtol += stepcoordtol

		N = len(fltrdhdu)

	fltrdcoordtbl = Table()
	fltrdcoordtbl.add_column(fltrdhdu['RAWX'], index = 0)
	fltrdcoordtbl.add_column(fltrdhdu['RAWY'], index = 1)
	fltrdcoordtbl.add_column(fltrdhdu['X'], index = 2)
	fltrdcoordtbl.add_column(fltrdhdu['Y'], index = 3)
	#print ( fltrdcoordtbl )

	#print ( (np.mean(fltrdhdu['RAWX']), np.mean(fltrdhdu['RAWY']) ), (np.mean(fltrdhdu['X']), np.mean(fltrdhdu['Y'])) )

	return [np.mean(fltrdhdu['X']), np.mean(fltrdhdu['Y'])]
    

def circ_reg_maker(regtype='circle', ra=500, dec=500, radius=18, outfile = "outfile.reg", a = None, b = None, angle=0.0) :
	if regtype == 'box' :
		flname = open(outfile, 'w')
		flname.write("# Region file format: DS9 version 4.1\n")
		text4regf = 'global color=green dashlist=8 3 width=1 font="helvetica'
		text4regf +=  '10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1'
		text4regf +=  ' delete=1 include=1 source=1\n'
        
		flname.write(text4regf)
		flname.write('image\n')
		flname.write('{}({}, {}, {}, {})'.format(regtype, ra, dec, a, b, angle))
		flname.close()

	elif regtype == 'circle' :
		flname = open(outfile, 'w')
		flname.write("# Region file format: DS9 version 4.1\n")
		text4regf = 'global color=green dashlist=8 3 width=1 font="helvetica'
		text4regf +=  '10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1'
		text4regf +=  ' delete=1 include=1 source=1\n'

		flname.write(text4regf)
		flname.write('image\n')
		flname.write('{}({}, {}, {})'.format(regtype, ra, dec, radius))
		flname.close()
	return outfile

def main() :

	import numpy as np
	import matplotlib.pyplot as plt
	import optparse, os, shutil, glob
	from astropy.io import fits
	from astropy.table import Table, Column, vstack
	import seaborn; seaborn.set()
	import scipy
	import scipy.signal
	import scipy.stats as st
	import pyregion, shutil

	usage = "usage: %prog [options] "
	parser = optparse.OptionParser(usage)

	parser.add_option("-e", "--evtflname", dest = "evtflname", 
						help = "The name of input events file \n Default" 
						+"is None which will end the procedure with no results", default = None)
						
	parser.add_option("-f", "--regfile", dest = "regfile", 
						help = "The name input region file \n Default" 
						+"is None which indicates full frame products", default = None)

	parser.add_option("-r", "--radreg", dest = "radreg", 
						help = "The string input to define the size of source region\n" 
						+ "The format is rad1nrad2nrad3..., where 'n' is separator. \n"
						+ "Write radn if only one radii is to be given \n"
						+ "defalut is '15' which reflects single radii of 15 arcmin", default = 15)

	parser.add_option("-a", "--autonbin", dest = "autoNbin", 
						help = "The nuber of bins for auto binning the curve", default = 30)

	parser.add_option("", "--outf", dest = "outf", help = "The Stem for output", default = None)
		
	parser.add_option("-u", "--chmin", dest = "chmin", 
						help = "Channel min cutoff\n Default Value = 31, i.e. 0.3 keV", default=31)
	parser.add_option("-v", "--chmax", dest = "chmax", 
						help = "Channel max cutoff\n Default Value = 701, i.e. 7.0 keV", default=701)

	parser.add_option("-g", "--grdflag", dest = "grdflag", 
						help = "Events Grade Selection \n Default Value = \"0-12\"", default="0-12")

	parser.add_option("-l", "--logname", dest = "logfile", 
						help = "The name of output logfile for debugging, \n"
						+ "Default Name = \"BkgAnaLogFile.log\"", default = "BkgAnaLogFile.log")

	parser.add_option("-o", "--fouttbl", dest = "FinalOutTable", 
						help = "The IPAC table name for output \n"
						+ "Default is 'None' which will append '.tbl'"
						+ "to input input events file name", default = None)

	parser.add_option("-t", "--regtype", dest = "regtype", 
						help = "The type of region selection i.e. circle or box \n"
						+ "Enter 'c' or 'b' for cicle or box \n Defalut is c (circle)", default = 'c')

	parser.add_option("", "--xyname", dest="xyname", 
						help = "The input for the XYNAME keyword for \"XSELECT\" "
						+ "to be used for products generation \n e.g., RAW or Sky \n"
						+ "Default Value for this parameter is \"sky\"", default = "sky")
						
	parser.add_option("", "--timefile", dest="timefile", 
						help = "The name of time file", default = None)
						
	parser.add_option("", "--intensity", dest="intensity", 
						help = "The intensity filter", default = None)	
										
	parser.add_option("", "--proddir", dest="proddir", 
						help = "The directory name to keep products", default = "Results")					

	(options, args) = parser.parse_args()

	evtflname  = options.evtflname
	outfstem =  options.outf
	grade_flag = options.grdflag
	regfile = options.regfile
	radreg = str(options.radreg)
	autoNbin = options.autoNbin
	timefile = options.timefile
	intensityFltr = options.intensity
	proddir = options.proddir
	
	chanmin, chanmax = int(options.chmin), int(options.chmax)
	outlogfile, FinalOutTable = options.logfile, options.FinalOutTable
	regtype, xyname = options.regtype, options.xyname
	
	inputpars = "\n---------Input Parameters----------------\n"
	for keywds in ["evtflname", "outfstem", "grade_flag", "regfile", "radreg",
					"autoNbin", "chanmin", "chanmax", "outlogfile", 
					"FinalOutTable", "regtype", "xyname", "timefile", "intensityFltr", 
					"proddir"]:
						
		inputpars += f"{keywds}: \t{eval(keywds)}\n"
	inputpars += "-------------------------\n"
	
	print(inputpars)	
	if xyname in ['raw', 'r', 'RAW', 'R']:
		xyname = "raw"
		imgflag = 'r'
	else:
		xyname = "sky"	
		imgflag = 's'
   
	curvebinsize = (float(autoNbin) * 2.3774729)/2 #TIMEDEL
    
	enstring = str(round(chanmin/100.,1)).split(".")[0]+"p"
	enstring += str(round(chanmin/100.,1)).split(".")[1]+"to"
	enstring += str(round(chanmax/100.,1)).split(".")[0]+"p"
	enstring += str(round(chanmax/100.,1)).split(".")[1]

	if evtflname == None :
		print ()
		evtflname = input("No events file is provided... Input Now!!")
		if not os.path.exists(evtflname):
			os._exit(199)

	if FinalOutTable == None :
		FinalOutTable = evtflname.split("_cl.evt")[0]+"out.tbl"

	if regtype == 'c' :
		regtype = 'circle'
	elif regtype == 'b' :
		regtype = 'box'
	else :
		print ("Enter proper regtype")
		os._exit()

	cwdir = os.getcwd()
	resdir = os.path.join(cwdir, proddir)
	if not os.path.exists(resdir) :
		os.makedirs(resdir)
		
	if outfstem == None:
		SrcNameOut =  evtflname.split("_cl.evt")[0]
	else:
		SrcNameOut = outfstem
		SrcNameOut = SrcNameOut.replace('.', 'p')
		
	SrcNameOut = f"{SrcNameOut}_{imgflag}_Gd{grade_flag}"
	
	if regfile != None:
		#running the product generator using the input region file...
		status = StdProductMaker(evtflname, regdir = "./", regfile = regfile, CHANMIN=chanmin, 
									CHANMAX=chanmax, productstem = SrcNameOut, 
									grade_flag = grade_flag, curvebinsize = curvebinsize, 
									xyname = xyname, intensity = intensityFltr,
									timefile = timefile)
		
		srcNameStem = os.path.splitext(SrcNameOut)[0]							
		oregfname = f"{srcNameStem}_o.reg"
		for ofiles in glob.glob(f"{srcNameStem}*"):
			if not os.path.splitext(ofiles)[-1] == ".evt":
				shutil.move(ofiles, os.path.join(resdir, ofiles))
			
		shutil.move(regfile, os.path.join(resdir, oregfname))							
		#os.system("mv {}* {}".format(SrcNameOut, resdir) )
									
	else:
		radvec = np.array(radreg.split("n"), float)
		for j in range(len(radvec)) :
			radii = radvec[j]  #in arcm
			radii_str = str(np.round(radii,2))
			radii_str = radii_str.replace(".","p")
			radii = round(radii*60.0/4.122, 1)  #in pixel

			if len(radii_str.split('.')) >1 :
				radii_str = radii_str.split('.')[0] + 'p' + radii_str.split('.')[1]

			if xyname.upper() == "RAW" :
				outregfile = circ_reg_maker(regtype = regtype, ra=300, dec=300, radius = radii, 
										outfile = "outfile.reg", a = None, b = None, angle=0.0)
										
				#print (f"Converted Centeral Position: {CenCoord}")						
			else :
				CenCoord = GetCentralcoordLookup(evtflname, rawx=300.0, rawy = 300.0)
				cenra, cendec = CenCoord
				print (f"Converted Centeral Position: {CenCoord}")
				outregfile = circ_reg_maker(regtype = regtype, ra=cenra, dec=cendec, radius = radii,
										outfile = "outfile.reg", a = None, b = None, angle=0.0)

			SrcNameOut2use = SrcNameOut.split("_Gd")[0]+ f"_Rd{radii_str}"
			SrcNameOut2use = f"{SrcNameOut2use}_Gd{grade_flag}"
			
			
			status = StdProductMaker(evtflname,regdir = "./", regfile = outregfile, CHANMIN=chanmin, 
										CHANMAX=chanmax, productstem = SrcNameOut2use, 
										grade_flag = grade_flag, curvebinsize = curvebinsize, 
										xyname = xyname, intensity = intensityFltr,
										timefile = timefile)
										
			srcNameStem = os.path.splitext(SrcNameOut2use)[0]
			
			#plotting complotter...
			ProfPlotterComb(srcNameStem, enstring, grade_flag, radii_str, figsize = [12, 4],
								resdir = proddir, imgflag = imgflag)
			print(inputpars)
			#lcfile
			lcfname = os.path.join(resdir, f"{SrcNameOut}.lc")

def ProfPlotterComb(srcNameStem, enstring, grade_flag, radii_str, figsize = [12, 4], 
					resdir = None, imgflag = 's'):
	import glob, os, shutil
	import matplotlib.pyplot as plt
	
	if resdir == None:
		resdir = "Results"
	if not os.path.exists(resdir):
		os.makedirs(resdir)
	
	pltdirs = os.path.join(resdir, "PlotDir")
	if not os.path.exists(pltdirs):
		os.makedirs(pltdirs)
	
	#pltdirs = os.path.join(resdir, pltdirs)	
						
	pltfiles = []
	windsze = [10,10]
	usefulfilelist = glob.glob(f"{srcNameStem}*")
	for ofiles in usefulfilelist:
		if not os.path.splitext(ofiles)[-1] == ".evt":
				
			if os.path.splitext(ofiles)[-1] == ".pha":
						
				phaof = f"{os.path.splitext(ofiles)[0]}_sp.png"
				s, ax3 = phaplotter(ofiles, grpmin = 30, outfile = phaof, figsize= windsze,
												phatyp = 'unbinned', bkgflag = False)
				pltfiles.append(phaof)
				#ax3add = ax3			
						
			if os.path.splitext(ofiles)[-1] == ".lc":
				lcof = f"{os.path.splitext(ofiles)[0]}_{enstring}_lc.png"
				x, ax2 = lcplotter(ofiles, binsize = 60, en = enstring, outfile = lcof,
										LCOutTable = os.path.splitext(lcof)[0]+".ipac",
										figsize = windsze)	
				pltfiles.append(lcof)		
				shutil.move(f"{os.path.splitext(lcof)[0]}.ipac", 
								os.path.join(pltdirs, f"{os.path.splitext(lcof)[0]}.ipac"))
															
			if os.path.splitext(ofiles)[-1] == ".img":
				imof = f"{os.path.splitext(ofiles)[0]}_{enstring}_im.png"
				x, ax1 = imgaeplotter(infile = ofiles, outfile = imof, figsize = windsze,
				                       imgflag = imgflag)
						
				pltfiles.append(imof)
											
			shutil.move(ofiles, os.path.join(resdir, ofiles))
									
	fig = plt.figure(figsize = figsize, constrained_layout=True)
	fig.tight_layout()
	#fig.add_subplot()
	#SORTING THE FILES in THE LIST
	pltfiles.sort(key=os.path.splitext)
				
	for i, s in enumerate(pltfiles):
		if len(pltfiles) > 3:
			fig.add_subplot(2, 2, i+1)
		else:
			fig.add_subplot(1, 3, i+1)
								
		plt.imshow(plt.imread(s))
		plt.axis('off')
		os.remove(s)
		
	combpltf = f"{srcNameStem}_{enstring}.png"	
	
	#fig.subplots_adjust(hspace=0.0, wspace=0.1)
	fig.savefig(combpltf, transparent=False, dpi=100,facecolor='w',edgecolor='w',
							orientation='portrait',	pad_inches=0.1)
							
	shutil.move(combpltf, os.path.join(pltdirs, combpltf))
	
	
def coordinateConverter(coordin = [15.5432, 15.329], coordtyp = 'deg2hhmmss'):

    import numpy as np
    import os, shutil
    print (f"You have entered coordinates : $\alpha$ = {coordin[0]} & $\dalta$ = {coordin[1]}")
    if coordtyp in ["DEG2hhmmss", "Deg2hhmmss", "DEg", "deg2hhmmss", "d2hhmmss", "d2hh", "d2h", "D2H", "D2h", "d2H"] :
        RAindeg = float(coordin[0])
        DECindeg = float(coordin[1])

        RAhh = int(RAindeg/15.)

        RAmmtmp = ( (RAindeg/15.)%1 )* 15.0 * 4.0

        RAmm = int(RAmmtmp)

        RAss = np.round((RAmmtmp%1) * 60.0, 2)

        DECdd = int(DECindeg)

        DECmmtmp =   60.0 * (abs(DECindeg)%1)

        DECmm = int(DECmmtmp) 

        DECss = np.round( (DECmmtmp%1) * 60., 2)

        RA2return = "{}:{}:{}".format(RAhh, RAmm, RAss)
        DEC2return = "{}:{}:{}".format(DECdd, DECmm, DECss)

    elif coordtyp in ["hhmmss2DEG", "hhmmss2Deg", "hhmmss2DEg", "hhmmss2deg", "hhmmss2d", "hh2deg", "h2d", "H2D", "HH2DD", "HHMMSS2DD"] :

        RAstr = str(coordin[0])
        DECstr = str(coordin[1])

        RAindeg = float(RAstr.split(":")[0]) * 15.0 + (float(RAstr.split(":")[1]) / 4.0) + (float(RAstr.split(":")[2]) / 240.)
        DECindeg = float(DECstr.split(":")[0]) + (float(DECstr.split(":")[1])/60.) + (float(DECstr.split(":")[2])/3600.)

        RA2return = RAindeg
        DEC2return = DECindeg

    else :
       print ("Enter proper tag for the format of convrsion")
       RA2return, DEC2return = None, None

    return [RA2return, DEC2return]
       
if __name__ == "__main__":
    main()
     

 
