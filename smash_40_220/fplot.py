#!/usr/bin/env python

'''

      Toy program to demo a query of the Data Lab TAP service and visualize
  the object density and associated CMD.  The user moves the cursor in the
  left-hand plot, the smaller plots track to show the catalog points zoomed
  in around the cursor, and the g-r vs g CMD for the points around the cursor,
  and a histogram of the differential values.

  A reticle in the zoomed position plot shows the size of the sampling cursor
  being used.  A mouseover in the CMD plot will show the corresponding point
  in the density plot with a highlight marker.  Holding the Shift key down 
  will suspend tracking to allow you to move the cursor to another window
  without changing the zoom/CMD display.

  A box may be dragged out in the CMD plot to select points corresponding to 
  specific color/magnitude region, the density plot will automatically be 
  updated.  Similarly, the contrast of the density plot can be adjusted by
  selecting a horizontal region of the histogram plot, e.g. to show only the
  values 3-sigma above the mean drag the cursor between zero (the difference
  from the mean) and 3 (3-sigma above the mean). 

  Command-line options:
       -h, --help              	# Help
       -d, --debug           	# Debug
       -v, --verbose         	# Verbose

       -c, --catalog		# Catalog to query
       -l, --load		# Load the FITS bintable
       -p, --ra			# Set RA position
       -r, --dec		# Set DEC position
       -s, --sz			# Set field size

       -S, --small_k		# Set small kernel size
       -B, --big_k		# Set big kernel size
       -F, --floor		# Set clipping floor (sigma)
       -C, --ceiling		# Set clipping ceiling (sigma)

  Additional work to pre-fetch/cache data, or use of different density tools
  will help speed and precision but is TBD.

  MJF  (2/15/16)

'''



import matplotlib.pyplot as plt
import numpy as np
from gavo import votable
from lxml import etree
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from matplotlib.widgets import Button, RadioButtons, CheckButtons
from matplotlib.widgets import Slider
from matplotlib.widgets import Cursor
from matplotlib.widgets import RectangleSelector
from matplotlib.widgets import SpanSelector
from astropy import convolution
from astropy.io import fits
from astropy import convolution
from astropy.table import Table
import matplotlib.patches as patches
import mpld3
import sys, time, math, getopt
import psycopg2                 # for database connectivity
import psycopg2.extras



################################################
# Some interesting initial positions
################################################

xpos   = 0
ypos   = 0

xsize  = 2.5		# initial field size in degrees
ysize  = 2.5


################################################
# Globals
################################################

ra0     = None		# Entire dataset
dec0    = None
data0   = None
g0      = None
gr0     = None
ra_idx  = None
dec_idx = None

fname   = None		# file name to load
field   = 169		# Hydra II

small_k = 40.0		# Differential kernel parameters
big_k   = 220.0
floor   = 2.0
ceiling = 999.0

use_tap	= 1

dmap    = None		# resultant density map

cursz  = 5.0		# mouseover cursor size (degrees)


# Data Lab TAP Service endpoint
accessURL = "http://zeus1.sdm.noao.edu:8080/ivoa-dal/tap"
        
# Direct database access connection
conn_str = "host='zeus1.sdm.noao.edu' dbname='tapdb' user='dlquery' password=''"
conn_str = "host='localhost' dbname='tapdb' user='dlquery' password=''"

# Default catalog
catalog = 'smash.blue_stars'


#  GETDATA -- Retrieve the data via TAP and put into a numpy array

def getData (field, catalog='smash.blue_stars'):
    global   xpos, ypos, xsize, ysize, use_tap

    #  Template query string based on (xpos,ypos) position.
    query = ('select ra_j2000,dec_j2000,gmag,rmag from '+ catalog + \
       ' where (fid = \'%d\' AND' \
       ' (gmag is not null) and ' + \
       ' (gmag between 10 and 23) and ((gmag-rmag) between -1.0 and 1.0))') % \
             field

    print "Getting data ....",
    sys.stdout.flush()

    start_time = time.time()

    if (use_tap == 1):
        raw = votable.ADQLSyncJob (accessURL, query).run().openResult()

        # We now need to parse the VOTable XML
        xml = etree.parse (raw)
        data = []
        for row in xml.findall('//TR'):
            data.append ([td.text for td in row.findall('TD')])
        data = np.array (data).T

        ra0 = data[0].astype('float64')
        dec0 = data[1].astype('float64')
        g0 = data[2].astype('float')
        r0 = data[3].astype('float')
        gr0 = np.subtract (g0,r0)
        print 'Got', data.shape[1], 'objects in %s seconds' % \
            (time.time()-start_time)

    else:
        try:
            conn = psycopg2.connect (conn_str)
            cursor = conn.cursor ()
        except psycopg2.Error as e:
            print 'Error: Cannot connect to database'
            sys.exit(2)

        cursor.execute (query)
        cursor.commit ()

        records = cursor.fetchall ()
        ra0 = records[1]
        dec0 = records[2]
        g0 = records[3]
        r0 = records[4]
        gr0 = np.subtract (g0,r0)

    xpos, ypos = np.mean (ra0), np.mean (dec0)
    xsize = (ra0.max() - ra0.min()) / 2.0
    ysize = (dec0.max() - dec0.min()) / 2.0

    return ra0, dec0, g0, gr0, data



#  READDATA -- Read the data from a static FITS BINTABLE (testing/demo mode)

def readData (fname):
    global   xpos, ypos, xsize, ysize
    global   ra0, dec0, g0, gr0

    start_time = time.time()
    print "Getting data in file '%s' ...." % fname,
    sys.stdout.flush()

    data = fits.getdata (fname)
    t = Table (data)

    ra0 = t['ra_j2000']
    dec0 = t['dec_j2000']
    #g0 = t['g']		# FIXME -- Need to account for null values
    #gr0 = t['g_r']
    g0 = t['gmag']		# FIXME -- Need to account for null values
    gr0 = t['gmag'] - t['rmag']

    print 'c/m  (%g,%g) (%g,%g)' % (g0.min(),g0.max(),gr0.min(),gr0.max())

    xpos, ypos = np.mean (ra0), np.mean (dec0)
    xsize = (ra0.max() - ra0.min()) / 2.0
    ysize = (dec0.max() - dec0.min()) / 2.0

    print 'Got', len(ra0), 'objects in %s seconds' % (time.time()-start_time)

    return ra0, dec0, g0, gr0, data




# DWARF_FILTER -- Ken Mighell's differential convolution filter.

def dwarf_filter (ra, dec, ra_index, dec_index, fwhm_small=2.0, fwhm_big=20):

    if (ra is None or dec is None):
        return

    x = (ra if ra_index is None else ra[ra_index])
    y = (dec if dec_index is None else dec[dec_index])


    start_time = time.time()  
    print "Computing differential convolution .... ",
    sys.stdout.flush()

    # Information about declination (y) [degrees]
    ymean = (y.min() + y.max()) / 2.0
    ydiff_arcmin = (y.max() - y.min()) * 60.0 # convert from degrees to arcmin

    # Information about right ascension (x) [degrees in time]:
    xdiff = x.max() - x.min() # angular separation [degrees (time)] 
    xmean = (x.min() + x.max())/2.0

    # convert from degrees in time to separation in angular degrees:
    xdiff_angular = (x.max() - x.min()) * np.cos (ymean*(np.pi/180.0))

    # convert from degress to arcmin
    xdiff_angular_arcmin = xdiff_angular * 60.0 

    # Get the number of one-arcmin pixels in the X and Y directions:
    nx = np.rint (xdiff_angular_arcmin).astype('int')
    ny = np.rint (ydiff_arcmin).astype('int')

    # Create a two-dimensional histogram of the raw counts:
    Counts, xedges, yedges  = np.histogram2d (x, y, (nx,ny) )
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    raw_hist = np.rot90 (Counts).copy() # hack around Pythonic weirdness

    # Make the small and big Gaussian kernels with a standard deviation
    # of the given FWHM in arcmin^2 pixels.
    kernel_small = convolution.Gaussian2DKernel (
        stddev = (fwhm_small / 2.35), factor=1)
    kernel_big = convolution.Gaussian2DKernel (
        stddev = (fwhm_big / 2.35), factor=1)

    # Compute the differential convolution kernels.
    conv_big = convolution.convolve (raw_hist, kernel_big)
    conv_small = convolution.convolve (raw_hist, kernel_small)
    conv_delta = conv_small - conv_big
    delta = conv_delta.copy()

    # Compute statistics and the floor
    mean = np.mean (delta, dtype='float64')
    sigma = np.std (delta, dtype='float64')
    median = np.median (delta)                       # not used
    floor = mean

    #print 'dwarf_filter: mean = %g  sigma = %g' % (mean, sigma)


    # Clip to specified limits.
    clipped = delta.copy()

    if dmap is not None:
        if dmap.lower == None:
            dmap.lower = mean + (2 * sigma)
        else:
            floor = dmap.lower
        clipped[ delta < dmap.lower ] = dmap.lower

        if dmap.upper != None:
            clipped[ delta > dmap.upper ] = dmap.upper
        else:
            dmap.upper = 999.0

        #print 'clip limits = (%g,%g)' % (dmap.lower, dmap.upper)
    else:
        clipped[ delta < floor ] = floor

    print " done in %s seconds\n" % (time.time() - start_time)

    # Return the computed fields.
    return raw_hist, extent, delta, clipped
pass




#  DENSMAP -- A Class to handle the density convolution map..

class densMap (object):
    """
        DENSPLOT -- A class to  ...
    """

    def __init__ (self, ra, dec, small_k, big_k):

        self.ra = ra
        self.dec = dec
        self.small_k = small_k
        self.big_k = big_k

        self.raw = None
        self.extent = None
        self.delta = None
        self.clipped = None

        self.lower = None
        self.upper = None

    def setData (self, ra, dec):
        self.ra = ra
        self.dec = dec

    def clipRefresh (self, lower, upper):
        self.clipped = dmap.delta.copy()

        sigma = np.std (self.clipped, dtype='float64')
        floor = (self.clipped.mean() if lower is None else lower)

        self.lower = lower * sigma

        self.clipped[ self.delta < self.lower ] = self.lower
        if upper is not None:
            self.upper = upper * sigma
            self.clipped[ self.delta > self.upper ] = self.upper

        print 'clipRefresh: sigma = %g  floor = %g  upper = %g' % \
            (sigma,self.lower,self.upper)

        densWindow.redraw() 			# Redraw the sub-plots

    def setHist (self, raw, extent, delta, clipped):
        self.raw = raw				# raw counts histogram
        self.extent = extent			# edges of the field
        self.delta = delta			# difference in convolved field
        self.clipped = clipped			# clipped map of difference



#  DENSPLOT -- A Class to handle the density window.

class densPlot (object):
    """
        DENSPLOT -- A class to  ...
    """

    def __init__ (self, ax):

        self.ax = ax
        self.canvas = ax.figure.canvas

        self.xdata = None
        self.ydata = None
        self.showContours = True
        self.showImage = True
        self.use_raw = True

        # Create a text box to display the density stats.
        at = AnchoredText ("Max=%8.5g Min=%8.5g" % (0.0, 0.0),
            prop=dict(size=8), frameon=True, loc=2)
        self.ax.add_artist (at)


    def setData (self, xdata, ydata):
        self.xdata = xdata
        self.ydata = ydata

    def setLabels (self, xlab, ylab):
        self.ax.set_xlabel (xlab)
        self.ax.set_ylabel (ylab)

    def drawContour (self, hist, extent):
        x = np.linspace (extent[0], extent[1], num=len(hist[0]))
        y = np.linspace (extent[3], extent[2], num=len(hist))
        CS = self.ax.contour (x, y, hist, 8, alpha=0.25)
        self.ax.clabel (CS, inline=1, fontsize=10)
        self.contours = CS
        self.canvas.draw()			# redraw the plot


    def update (self):
        #global raw, extent, delta, clipped

        # Create the main density display.
        raw, extent, delta, clipped = dwarf_filter (ra0, dec0,
            ra_idx, dec_idx, fwhm_small=small_k, fwhm_big=big_k)
        dmap.setHist (raw, extent, delta, clipped)
        self.redraw()

    def redraw (self):
        self.ax.clear()

        self.ax.set_xlim (xpos + xsize, xpos - xsize)
        self.ax.set_ylim (ypos - ysize, ypos + ysize)

        at = AnchoredText ("Max=%8.5g Min=%8.5g  Npts=%d" % \
            (dmap.clipped.max(), dmap.clipped.min(),len(ra0[ra_idx])),
            prop=dict(size=8), frameon=True, loc=2)
        self.ax.add_artist (at)

        if self.showImage:
            self.ax.imshow (dmap.clipped, extent=dmap.extent, 
                interpolation='nearest', cmap='Greys')
        if self.showContours:
            if self.use_raw:
                self.drawContour (dmap.raw, dmap.extent)
            else:
                self.drawContour (dmap.clipped, dmap.extent)
        self.canvas.draw () 			# Redraw the sub-plots



#  CMDPLOT -- A Class to handle the CMD window.

class cmdPlot (object):
    """
        CMDPLOT -- A class to  ...
    """

    def __init__ (self, ax):
        global gr0, g0

        self.ax = ax
        self.canvas = ax.figure.canvas

        self.color   = None
        self.mag     = None

        self.ax.set_xlim (-0.5, 1.0)
        #self.ax.set_ylim (g0.max()+0.25, g0.min()-0.25)
        self.ax.set_ylim (min(25.0,g0.max()+0.25), 14.0)

    def setData (self, color, mag):
        self.color = color
        self.mag = mag

        self.ax.set_xlim (-0.5, 1.0)
        self.ax.set_ylim (min(25.0,mag.max()+0.25), 14.0)

    def setLabels (self, xlabel, ylabel):
        self.ax.set_xlabel (xlabel)
        self.ax.set_ylabel (ylabel)

    def update (self, x, y):
        global ra, dec, cursz


        try:
            # Find all points within some specified radius
            dist = np.hypot (x - ra0, y - dec0)
            idx = dist < cursz
            gr, g = gr0[idx], g0[idx]

            # Redraw the CMD plot
            self.ax.clear()
            self.ax.scatter (gr, g, s=1, marker='.', cmap='grey', alpha=0.09)
            at = AnchoredText ("Npts = %d" % len(gr), prop=dict(size=8), 
                frameon=True, loc=2)
            self.ax.add_artist (at)
           
            self.ax.set_xlim (-0.5, 1.0)
            self.ax.set_ylim (min(25.0,self.mag.max()+0.25), 14.0)

            self.canvas.draw () 		# Redraw the sub-plots

        except:
            # Just ignore any exceptions for now
            pass





############################
# Task main()
############################

if __name__ == '__main__':

    # Process commandline args.
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hdvtl:c:f:r:d:s:S:B:F:C:",
            ["help","debug","verbose","load","catalog=", "field=",
	     "ra=","dec=","sz=", "small_k=", "big_k=", "floor=", "ceiling="])

    except getopt.GetoptError:
        print 'Error processing arguments:'
        Usage ()
        sys.exit (1)


    for opt, arg in opts:
       if opt in ("-h", "--help"):              # Help
           Usage ()
           sys.exit ()
       elif opt in ("-d", "--debug"):           # Debug
           debug = 1
       elif opt in ("-v", "--verbose"):         # Verbose
           verbose = 1
       elif opt in ("-t", "--tap"):             # TAP interface
           use_tap = 1

       elif opt in ("-c", "--catalog"):		# Catalog to query
           catalog = arg
       elif opt in ("-f", "--field"):		# SMASH Field Number
           field = int (arg)
       elif opt in ("-l", "--load"):		# Load the FITS bintable
           fname = arg
       elif opt in ("-p", "--ra"):		# Set RA position
           xpos = float (arg)
       elif opt in ("-r", "--dec"):		# Set DEC position
           ypos = float (arg)
       elif opt in ("-s", "--sz"):		# Set field size
           xsize = ysize = float (arg)

       elif opt in ("-S", "--small_k"):		# Set small kernel size
           small_k  = float (arg)
       elif opt in ("-B", "--big_k"):		# Set big kernel size
           big_k  = float (arg)
       elif opt in ("-F", "--floor"):		# Set clipping floor (sigma)
           floor  = float (arg)
       elif opt in ("-C", "--ceiling"):		# Set clipping ceiling (sigma)
           ceiling  = float (arg)



    # Get the data and compute the differential convolution density
    # histograms.
    if fname is not None:
        ra0, dec0, g0, gr0, data0 = readData ( fname )

    elif xpos is not None and ypos is not None:
        ra0, dec0, g0, gr0, data0 = getData ( field, catalog )

    ra_idx = ra0 < 361				# Initialize for all values
    dec_idx = dec0 < 91


    # Initialize the density map.  [FIXME:  need to turn this into object]
    dmap = densMap (ra0, dec0, small_k, big_k)

    raw, extent, delta, clipped = dwarf_filter (dmap.ra, dmap.dec,
            ra_idx, dec_idx, fwhm_small=dmap.small_k, fwhm_big=dmap.big_k)
    dmap.setData (ra0, dec0)
    dmap.setHist (raw, extent, delta, clipped)


    # Create the subplot windows.
    fig, (axcmd, axdens) = plt.subplots (1, 2, figsize=(12,5))

    # Create the density convolution subplot window.
    densWindow = densPlot (axdens)
    densWindow.setLabels ('RA (deg)', 'Dec (deg)')
    densWindow.setData (ra0, dec0)
    densWindow.redraw ()
    axcmd.set_aspect (1.0)

    # Create a plot of the CMD centered on the initial position.
    cmdWindow = cmdPlot (axcmd)
    cmdWindow.setLabels ("(g-r)", "g")
    cmdWindow.setData (gr0, g0)
    cmdWindow.update (np.mean(ra0), np.mean(dec0))

    ratio_default=(axdens.get_xlim()[1]-axdens.get_xlim()[0]) / (axdens.get_ylim()[1]-axdens.get_ylim()[0])
    axdens.set_aspect (abs(ratio_default))
    ratio_default=(axcmd.get_xlim()[1]-axcmd.get_xlim()[0]) / (axcmd.get_ylim()[1]-axcmd.get_ylim()[0])
    axcmd.set_aspect (abs(ratio_default))
    fig.canvas.draw ()

    # Finalize the plot and display
    #plt.tight_layout (pad=2.0)
    plt.suptitle ("Field %d" % field, fontsize=20)
    #plt.show ()
    plt.savefig ("Field%d.jpeg" % field)

