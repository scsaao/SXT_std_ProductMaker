# SXT_std_ProductMaker
A handy command line tool create standard science products (image, lightcurves and spectra)

### Pre-Requisites: 
              #Python (>=3.5)
              #astropy
              numpy
              matplotlib
              shutil
              seaborn
              optparse
              scipy
              glob
              pyregion

### How To run
python src_bkg/SXT_Std_prod_Maker_v02.py -h <br />
Usage: SXT_Std_prod_Maker_v02.py [options] 

Options:
  * -h, --help            show this help message and exit <br />
  * -e EVTFLNAME, --evtflname=EVTFLNAME <br />
  The name of input events file, Default is <b> None</b> which will end the procedure with no results <br />
  * -f REGFILE, --regfile=REGFILE <br />
                        The name input region file   Defaultis None which indicates full frame products <br />
  * -r RADREG, --radreg=RADREG <br />
                        The string input to define the size of source region
                        The format is rad1nrad2nrad3..., where 'n' is
                        separator.  Write radn if only one radii is to be
                        given  defalut is '15' which reflects single radii of
                        15 arcmin <br />
  * -a AUTONBIN, --autonbin=AUTONBIN <br />
                        The nuber of bins for auto binning the curve <br />
  * --outf=OUTF      <br />     The Stem for output
  * -u CHMIN, --chmin=CHMIN <br />
                        Channel min cutoff  Default Value = 31, i.e. 0.3 keV 
  * -v CHMAX, --chmax=CHMAX <br />
                        Channel max cutoff  Default Value = 701, i.e. 7.0 keV
  * -g GRDFLAG, --grdflag=GRDFLAG <br />
                        Events Grade Selection   Default Value = "0-12"
  * -l LOGFILE, --logname=LOGFILE <br />
                        The name of output logfile for debugging,  Default
                        Name = "BkgAnaLogFile.log"
  * -o FINALOUTTABLE, --fouttbl=FINALOUTTABLE <br />
                        The IPAC table name for output  Default is 'None'
                        which will append '.tbl'to input input events file
                        name
  * -t REGTYPE, --regtype=REGTYPE <br />
                        The type of region selection i.e. circle or box  Enter
                        'c' or 'b' for cicle or box   Defalut is c (circle)
  * --xyname=XYNAME    <br />   The input for the XYNAME keyword for "XSELECT" to be
                        used for products generation   e.g., RAW or Sky
                        Default Value for this parameter is "sky"
  * --timefile=TIMEFILE   <br /> The name of time file
  * --intensity=INTENSITY <br />
                        The intensity filter
  * --proddir=PRODDIR   <br />  The directory name to keep products
