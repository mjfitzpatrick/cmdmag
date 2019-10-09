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

       -R, --reticle		# Set cursor size
       -Z, --zoomsz		# Set size of zoom box


  Density-plot Commands:
	+		    	# zoom in on field
	-		    	# zoom out on field
	h / left-arrow	    	# pan left (get new data)
	j / down-arrow	    	# pan down (get new data)
	k / up-arrow	    	# pan up (get new data)
	l / right-arrow	    	# pan right (get new data)

        Shift		    	# freeze motion tracking in density plot
        Shift-release	    	# enable motion tracking in density plot

  Zoom-plot Commands:
	< / >		    	# change sample cursor by a factor of 2
	, / .		    	# change sample cursor by a factor of 1.5
	+		    	# increase zoom factor
	-		    	# decrease zoom factor


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


################################################
# Some interesting initial positions
################################################

xpos   = 202.75655
ypos   = -8.09553

xpos   = 34.81666667	# Segue-2
ypos   = 20.17527777

xpos   = 168.3716667	# Leo - 2
ypos   = 22.15472222

xpos   = 263.5 		# DECaLS DR2 Field
ypos   = 26.9

xpos   = 229.015	# Pal-5
ypos   = -0.093

xpos   = 260.31		# DECaLS DR2 Field, possible horiz. branch?
ypos   = 25.46

xpos   = 247.758333 	# Hercules dwarf
ypos   = 12.775

xpos   = 185.4254	# Hydra-II
ypos   = -31.98523


xsize  = 0.5		# initial field size in degrees
ysize  = 0.5

xsize  = 0.75		# initial field size in degrees
ysize  = 0.75

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

small_k = 2.0		# Differential kernel parameters
big_k   = 20.0
floor   = 0.0
ceiling = 2.5

dmap    = None		# resultant density map

cursz  = 0.1		# mouseover cursor size (degrees)
zoomsz = 0.15		# initi size of zoom box (degrees)


# Data Lab TAP Service endpoint
accessURL = "http://dldb1.sdm.noao.edu:8080/ivoa-dal/tap"
accessURL = "http://zeus1.sdm.noao.edu/tap"

# Default to use the LS DR2 catalog
catalog = 'lsdr2.star'
catalog = 'smash.blue_stars'




#  GETDATA -- Retrieve the data via TAP and put into a numpy array

def getData (catalog='lsdr2.star'):

    print "catalog is " + catalog

    #  Template query string based on (xpos,ypos) position.
    if catalog == 'lsdr2.star':
        query = ('select ra_j2000,dec_j2000,g,g_r,c30star from '+ catalog + \
            ' where ' + \
            ' (ra_j2000 between %g and %g) and ' + \
            ' (dec_j2000 between %g and %g) and ' + \
            ' (g is not null) and ' + \
            #' (c30star < 40) and ' + \
            #' (g between 20.5 and 23.5) and (g_r between 0.25 and 0.5)') % \
            ' (g between 1 and 23) and (g_r between -0.5 and 2.0)') % \
                (xpos-xsize, xpos+xsize, ypos-ysize, ypos+ysize)
    else:
        query = ('select ra_j2000,dec_j2000,gmag,(gmag-rmag) as g_r,rmag,sharp '+
            ' from '+ catalog + \
            ' where ' + \
            ' (ra_j2000 between %g and %g) and ' + \
            ' (dec_j2000 between %g and %g) and ' + \
            ' (gmag is not null) and (rmag is not null) and ' + \
            ' (gmag between 10 and 23) and ' + \
            ' ((gmag-rmag) between -0.5 and 1.0)') % \
                (xpos-xsize, xpos+xsize, ypos-ysize, ypos+ysize)

    print "query is " + query
    print "Getting data ....",
    sys.stdout.flush()

    start_time = time.time()
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
    gr0 = data[3].astype('float')

    print 'Got', data.shape[1], 'objects in %s seconds' % \
        (time.time()-start_time)

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
    #g0 = t['g']			# FIXME -- Need to account for null values
    #gr0 = t['g_r']
    g0 = t['gmag']			# FIXME -- Need to account for null values
    gr0 = t['gmag'] - t['rmag']

    print type(ra0)
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
            dmap.lower = mean
        else:
            floor = dmap.lower
        clipped[ delta < dmap.lower ] = dmap.lower

        if dmap.upper != None:
            clipped[ delta > dmap.upper ] = dmap.upper
        else:
            dmap.upper = 999.0

        print 'clip limits = (%g,%g)' % (dmap.lower, dmap.upper)
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
        self.use_raw = False

        self.track_cursor = True

        # Create a text box to display the density stats.
        at = AnchoredText ("Max=%8.5g Min=%8.5g" % (0.0, 0.0),
            prop=dict(size=8), frameon=True, loc=2)
        self.ax.add_artist (at)

        # Create a marker to map a point in the CMD to it's place on the
        # plot.
        self.marker = patches.Ellipse ((xpos, ypos), width=0.05, height=0.05,
            angle=0, alpha=0.9, facecolor="red")
        self.ax.add_patch (self.marker)

        # Attach the event handlers.
        self.canvas.mpl_connect ('motion_notify_event', self.motionCallback)
        self.canvas.mpl_connect ('key_press_event', self.onpressCallback)
        self.canvas.mpl_connect ('key_release_event', self.onreleaseCallback)
        self.canvas.mpl_connect ('axes_leave_event', self.onleaveCallback)


    def onpressCallback (self, event):
        global xpos, ypos, xsize, ysize, cursz, zoomsz
        global ra0, dec0, g0, gr0, data0, ra_idx, dec_idx
        global small_k, big_k, ra_idx, dec_idx


        if (event.inaxes is not self.ax):
            return

        key = event.key
        if (key == 'up' or key == 'k'):		# Pan up
            ypos += 1.5 * ysize
        elif (key == 'down' or key == 'j'):	# Pan down
            ypos -= 1.5 * ysize
        elif (key == 'left' or key == 'h'):	# Pan left
            xpos -= 1.5 * xsize
        elif (key == 'right' or key == 'l'):	# Pan right
            xpos += 1.5 * xsize
        elif (key == 'c'):			# Center on current pos
            xpos, ypos = event.xdata, event.ydata

        elif (key == '('):			# Decrease fwhm_big
            big_k = max (5.0*small_k, big_k - 5)
        elif (key == ')'):			# Increase fwhm_big
            big_k += 5
        elif (key == '['):			# Decrease fwhm_small
            small_k = max (1.0, small_k - 2)
        elif (key == ']'):			# Increase fwhm_small
            small_k += 1

        elif (key == '<'):                      # Decrease cursor size * 2
            zoomWindow.setCursorSize (zoomWindow.cursz / 2.0)
        elif (key == '>'):                      # Increase cursor size * 2
            zoomWindow.setCursorSize (zoomWindow.cursz * 2.0)
        elif (key == ','):                      # Decrease cursor size * 1.5
            zoomWindow.setCursorSize (zoomWindow.cursz / 1.5)
        elif (key == '.'):                      # Increase cursor size * 1.5
            zoomWindow.setCursorSize (zoomWindow.cursz * 1.5)
        elif (key == '+'):                      # Decrease cursor size * 1.5
            zoomWindow.setZoomSize (zoomWindow.zoomsz * 1.5)
        elif (key == '-'):                      # Increase cursor size * 1.5
            zoomWindow.setZoomSize (max (zoomWindow.zoomsz / 1.5, 0.15))

        elif (key == 'shift'):			# Suspend cursor tracking
            self.track_cursor = False
            return
        else:
            return

        print 'Pos (%g,%g) : Size (%g,%g) : Cursor = %g : Zoom = %g' % \
            (xpos,ypos,xsize,ysize,cursz,zoomsz)

        if (key in ('h','j','k','l','left','down','up','right','-','+')):
            # Get the data for the new field.
            ra0, dec0, g0, gr0, data0 = getData ( catalog )

            ra_idx = ra0 < 361			# reset working dataset
            dec_idx = dec0 < 91
            self.update ()


        # Refresh the main density plot.
        raw, extent, delta, clipped = dwarf_filter (ra0, dec0, 
            ra_idx, dec_idx, fwhm_small=small_k, fwhm_big=big_k)
        dmap.setData (ra0[ra_idx], dec0[dec_idx])
        dmap.setHist (raw, extent, delta, clipped)

        # Update the zoomed data.
        zoomWindow.setData (ra0[ra_idx], dec0[dec_idx])
        zoomWindow.draw ()

        fig.canvas.draw () 			# Redraw the sub-plots


    def onreleaseCallback (self, event):
        if (event.inaxes is not self.ax):
            return
        if (event.key == 'shift'):		# resume motion tracking
            self.track_cursor = True

    def onleaveCallback (self, event):
        if (event.inaxes is not self.ax):
            return
        self.track_cursor = True		# resume motion tracking

    def motionCallback (self, event):
        if (event.inaxes is not self.ax):
            return
        if (self.track_cursor == False):
            return

        try:
            # Redraw the CMD and zoom plots
            cmdWindow.update (event.xdata, event.ydata) 		
            zoomWindow.setPos (event.xdata, event.ydata)
        except:
            # Just ignore any exceptions for now
            pass


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
        self.marker  = None

        self.ax.set_xlim (-0.5, 2.0)
        #self.ax.set_ylim (g0.max()+0.25, g0.min()-0.25)
        self.ax.set_ylim (max(25.0,g0.max()+0.25), 10.0)

        self.cmd_select = RectangleSelector (self.ax, 
                                           self.cmd_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # don't use middle
                                           minspanx=5, minspany=5,
                                           spancoords='pixels')

        # Attach the event handlers.
        self.canvas.mpl_connect ('motion_notify_event', self.motionCallback)


    def cmd_select_callback (self, eclick, erelease):
        global ra0, dec0, ra_idx, dec_idx

        # eclick and erelease are the press and release events
        x1, y1 = eclick.xdata, eclick.ydata
        x2, y2 = erelease.xdata, erelease.ydata

        print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))

        # Add a marker to indicate the selected region
        if (self.marker is not None):
            self.marker.set_visible (False)
        self.marker = patches.Rectangle ((x1, y1), 
                                         abs(x2-x1), abs(y2-y1), 
                                         color='red', alpha=0.25)
        self.ax.add_patch (self.marker)
        self.canvas.draw ()

        idx_x = np.logical_and (gr0 > x1, gr0 < x2)
        idx_y = np.logical_and (g0 > y1, g0 < y2)
        data_idx = np.logical_and (idx_x, idx_y)

        ra_idx = data_idx
        dec_idx = data_idx

        raw, extent, delta, clipped = dwarf_filter (ra0, dec0, 
            ra_idx, dec_idx, fwhm_small=small_k, fwhm_big=big_k)
        dmap.setData (ra0[ra_idx], dec0[dec_idx])
        dmap.setHist (raw, extent, delta, clipped)
        densWindow.redraw()


    def keypressCallback (self, event):
        if (event.inaxes is not self.ax):
            return

        # Not yet implemented
        pass


    def motionCallback (self, event):
        global ra0, dec0, marker

        if (event.inaxes is not self.ax):
            return

        try:
            # Find all points within some specified radius
            x, y = event.xdata, event.ydata
            dist = np.hypot (x - self.color, y - self.mag)
            indx = dist.argmin()
            dx, dy = ra[indx], dec[indx]

            # Redraw the density marker at the location corresponding to
            # the nearest point in the CMD.
            densWindow.marker.center = dx, dy
            self.canvas.draw ()
            fig.canvas.draw ()

        except:
            # Just ignore any exceptions for now
            pass


    def setData (self, color, mag):
        self.color = color
        self.mag = mag

        self.ax.set_xlim (-0.5, 2.0)
        self.ax.set_ylim (max(25.0,mag.max()+0.25), 10.0)


    def setLabels (self, xlabel, ylabel):
        self.ax.set_xlabel (xlabel)
        self.ax.set_ylabel (ylabel)


    def update (self, x, y):
        global ra, dec, cursz


        if (densWindow.track_cursor == False):
            return

        try:
            # Find all points within some specified radius
            dist = np.hypot (x - ra0, y - dec0)
            idx = dist < cursz
            gr, g = gr0[idx], g0[idx]

            # Redraw the CMD plot
            self.ax.clear()
            self.ax.scatter (gr, g, s=1, marker='.', cmap='jet')
            at = AnchoredText ("Npts = %d" % len(gr), prop=dict(size=8), 
                frameon=True, loc=2)
            self.ax.add_artist (at)
           
            self.ax.set_xlim (-0.5, 2.0)
            self.ax.set_ylim (max(25.0,self.mag.max()+0.25), 10.0)

            densWindow.marker.center = x, y

            self.canvas.draw () 		# Redraw the sub-plots
        except:
            # Just ignore any exceptions for now
            pass





#  ZOOMPLOT -- A Class to handle the position zoom window.

class zoomPlot (object):
    """
        ZOOMPLOT -- A class to display a magnifier window zoomed in on a
        region around the cursor in another plot window.

        Key Bindings:
              <		Decrease cursor size by factor of 2
              >		Increase cursor size by factor of 2
              ,		Decrease cursor size by factor of 1.5
              .		Increase cursor size by factor of 1.5
              +		Increase zoom area
              -		Decrease zoom area
    """

    def __init__ (self, ax):
        global xpos, ypos, cursz, zoomsz

        self.ax = ax
        self.canvas = ax.figure.canvas

        self.cursz   = cursz
        self.zoomsz  = zoomsz
        self.xdata   = None
        self.ydata   = None
        self.cur_x   = xpos
        self.cur_y   = ypos


        # Create a reticle overlay to show the sampling cursor size.
        self.reticle = patches.Ellipse ((xpos, ypos), 
            width=(2 * cursz), height=(2.66 * cursz),
            angle=0, alpha=0.1, facecolor="cyan")

        # Attach the event handlers.
        self.canvas.mpl_connect ('key_press_event', self.keypressCallback)


    def keypressCallback (self, event):
        global  cursz

        if (event.inaxes is not self.ax):
            return

        key = event.key
        if (key == '<'):                        # Decrease cursor size * 2
            self.setCursorSize (self.cursz / 2.0)
        elif (key == '>'):                      # Increase cursor size * 2
            self.setCursorSize (self.cursz * 2.0)
        elif (key == ','):                      # Decrease cursor size * 1.5
            self.setCursorSize (self.cursz / 1.5)
        elif (key == '.'):                      # Increase cursor size * 1.5
            self.setCursorSize (self.cursz * 1.5)

        elif (key == '+'):                      # Zoom In
            self.setZoomSize (self.zoomsz * 1.5)
        elif (key == '-'):                      # Zoom Out
            self.setZoomSize (max (self.zoomsz / 1.5, 0.15))

        elif (key == 'shift'):
            return
 

    def setData (self, x, y):
        self.xdata = x
        self.ydata = y


    def setCursorSize (self, size):
        global  cursz, zoomsz

        size = max (size, 0.01)			# enforce min cursor size
        cursz = size
        self.cursz = size

        # Note: The reticle size is adjusted for the aspect ratio of the plot.
        self.reticle.width = (2 * size)
        self.reticle.height = (2.66 * size)
        self.reticle.center = self.cur_x, self.cur_y # update reticle position

        zoomsz = max (0.15, cursz * 1.5)
        self.zoomsz = max (0.15, cursz * 1.5)

        self.ax.set_xlim (self.cur_x + zoomsz, self.cur_x - zoomsz)
        self.ax.set_ylim (self.cur_y - zoomsz, self.cur_y + zoomsz)

        cmdWindow.update(self.cur_x, self.cur_y)

        self.canvas.draw()			# redraw the plot


    def setZoomSize (self, size):
        global xpos, ypos, cursz, zoomsz

        zoomsz = size
        self.zoomsz = size

        # Flip RA for std orientation
        self.ax.set_xlim (self.cur_x + zoomsz, self.cur_x - zoomsz) 
        self.ax.set_ylim (self.cur_y - zoomsz, self.cur_y + zoomsz)

        self.canvas.draw()			# redraw the plot


    def setPos (self, x, y):
        self.cur_x = x
        self.cur_y = y

        # Redraw the zoom plot (flip RA for std orientation)
        self.ax.set_xlim (self.cur_x + zoomsz, self.cur_x - zoomsz)
        self.ax.set_ylim (self.cur_y - zoomsz, self.cur_y + zoomsz)
        self.reticle.center = x, y		# update reticle position
        self.canvas.draw()			# redraw the plot


    def setLabels (self, xlab, ylab):
        self.ax.set_xlabel (xlab)
        self.ax.set_ylabel (ylab)


    def draw (self):
        global xpos, ypos, cursz

        if (self.xdata is None) or (self.ydata is None):
            return

        # Create the zoomed positional plot.  Note the reticle is adjusted for
        # the aspect ratio of the subplot.
        self.ax.scatter (self.xdata, self.ydata, s=9, alpha=0.6, 
            marker='.', cmap='Greys')

        print "pos=(%g,%g) cursz=%g zoom=%g" % (xpos, ypos, cursz, zoomsz)

        self.ax.set_xlim (self.cur_x + self.zoomsz, self.cur_x - self.zoomsz) 
        self.ax.set_ylim (self.cur_y - self.zoomsz, self.cur_y + self.zoomsz)

        self.reticle = patches.Ellipse ((self.cur_x, self.cur_y), 
            width=(2 * cursz), height=(2.66 * cursz),
            angle=0, alpha=0.1, facecolor="cyan")
        self.ax.add_patch (self.reticle)

        atz = AnchoredText ("Zoom=%g Cursor=%g" % (zoomsz, cursz),
            prop=dict(size=8), frameon=True, loc=2)
        self.ax.add_artist (atz)


    def drawContour (self, hist, extent):
        x = np.linspace (extent[0], extent[1], num=len(hist[0]))
        y = np.linspace (extent[3], extent[2], num=len(hist))
        CS = self.ax.contour (x, y, hist)
        self.ax.clabel (CS, inline=1, fontsize=10)
        self.canvas.draw()			# redraw the plot




#  HISTPLOT -- A Class to handle the position zoom window.

class histPlot (object):
    """
        HISTPLOT -- A class to display a histogram window.

        Key Bindings:
    """

    def __init__ (self, ax, data):

        self.ax = ax
        self.canvas = ax.figure.canvas

        self.ax.get_yaxis().set_visible (False)

        self.marker  = None
        self.lower   = None
        self.upper   = None
        self.data    = data

        dv = np.ndarray.flatten (data)
        self.mean = np.mean (dv,dtype='float64')
        self.sigma = np.std (dv,dtype='float64')
        self.median = np.median (dv)

        self.lower = self.mean
        self.upper = self.mean + (2 * self.sigma)

        # Draw the histogram
        sighist = (dv-self.mean) / self.sigma
        self.ax.hist (sighist, 20, log=True, color='lightgrey')

        # FIXME
        self.ax.set_xlabel ('Sigma')
        self.ax.set_xlim (max(-6.0,data.min()), min(10.0,data.max()))

        self.ax.axvline (x=self.mean, color='red')
        self.ax.axvline (x=-3.0, color='cyan')
        self.ax.axvline (x=3.0, color='cyan')

        # set useblit True on gtkagg for enhanced performance
        self.span = SpanSelector (self.ax, self.onselect, 
                         'horizontal', useblit=True,
                         rectprops=dict(alpha=0.5, facecolor='yellow'))

        self.canvas.draw()


    def onselect (self, xmin, xmax):
        self.lower = max(self.mean,xmin)
        self.lower = xmin
        self.upper = xmax

        # Add a marker to indicate the selected region
        if (self.marker is not None):
            self.marker.set_visible (False)
        self.marker = patches.Rectangle ((self.lower, 0.0), 
                                         abs(xmax-self.lower), 99999.0, 
                                         color='red', alpha=0.25)
        self.ax.add_patch (self.marker)

        dmap.clipRefresh (xmin, xmax)




#  CONTROLFRAME -- A class to manage the control widgets.

class controlFrame (object):
    """
        CONTROLFRAME -- A class to manage the control widgets.
    """

    def __init__ (self):

        self.smallS = Slider (plt.axes ([0.1, 0.11, 0.325, 0.03]),
            'Small-K', 1, 40, valinit=small_k)
        self.bigS = Slider (plt.axes ([0.1, 0.075, 0.325, 0.03]),
            'Big-K', 1, 99, valinit=big_k)
        self.countS = Slider (plt.axes ([0.1, 0.04, 0.16, 0.03]),
            'C30 Stars', 10, 99, valinit=40)

        self.reset_button = Button (plt.axes([0.315, 0.04, 0.05, 0.03]),
            'Reset', color='white', hovercolor='yellow')
        self.apply_button = Button (plt.axes([0.375, 0.04, 0.05, 0.03]),
            'Apply', color='white', hovercolor='yellow')

        self.up_button = Button (plt.axes([0.25, 0.925, 0.05, 0.03]),
              'N', color='white', hovercolor='yellow')
        self.down_button = Button (plt.axes([0.25, 0.855, 0.05, 0.03]),
              'S', color='white', hovercolor='yellow')
        self.left_button = Button (plt.axes([0.22, 0.89, 0.05, 0.03]),
              'E', color='white', hovercolor='yellow')
        self.right_button = Button (plt.axes([0.28, 0.89, 0.05, 0.03]),
              'W', color='white', hovercolor='yellow')

        self.zoom_in = Button (plt.axes([0.375, 0.925, 0.095, 0.03]),
              'Zoom In', color='white', hovercolor='yellow')
        self.zoom_out = Button (plt.axes([0.375, 0.89, 0.095, 0.03]),
              'Zoom Out', color='white', hovercolor='yellow')
        self.center = Button (plt.axes([0.375, 0.855, 0.095, 0.03]),
            'Center', color='white', hovercolor='yellow')

        self.do_contour = CheckButtons (plt.axes([0.055, 0.855, 0.115, 0.095]),
              ('Contours','Image','Raw Bkg'), (True,True,False))

        # Bind the slider control methods.
        self.smallS.on_changed (self.smallk_update)
        self.bigS.on_changed (self.bigk_update)
        self.countS.on_changed (self.counts_update)

        # Bind the Pan/Zoom control methods.
        self.zoom_in.on_clicked (self.zoomIn)
        self.zoom_out.on_clicked (self.zoomOut)
        self.center.on_clicked (self.center)
        self.apply_button.on_clicked (self.apply)
        self.reset_button.on_clicked (self.reset)

        self.up_button.on_clicked (self.panUp)
        self.down_button.on_clicked (self.panDown)
        self.left_button.on_clicked (self.panLeft)
        self.right_button.on_clicked (self.panRight)

        # Bind the density plot options methods.
        self.do_contour.on_clicked (self.plotOpts)


    def panUp (self, event):
        global ypos, ysize
        ypos += 1.5 * ysize
        densWindow.update ()

    def panDown (self, event):
        global ypos, ysize
        ypos -= 1.5 * ysize
        densWindow.update ()

    def panLeft (self, event):
        global xpos, xsize
        xpos += 1.5 * xsize
        densWindow.update ()

    def panRight (self, event):
        global xpos, xsize
        xpos -= 1.5 * xsize
        densWindow.update ()

    def zoomIn (self, event):
        global xsize, ysize
        xsize /= 2.0
        ysize /= 2.0
        densWindow.update ()

    def zoomOut (self, event):
        global xsize, ysize
        xsize *= 2.0
        ysize *= 2.0
        densWindow.update ()

    def smallk_update (self, val):
        global small_k
        small_k = val

    def bigk_update (self, val):
        global big_k
        big_k = val

    def counts_update (self, val):
        print "counts val = ",		# NYI
        print val
        pass

    def apply (self, event):
        densWindow.update ()

    def center (self, event):
        pass				# NYI

    def reset (self, event):
        # ....reset cursor/zoom size to defaults
        self.smallS.reset()		# reset sliders
        self.bigS.reset()
        self.countS.reset()
        pass

    def plotOpts (self, label):
        if label == 'Contours':
            densWindow.showContours = (not densWindow.showContours)
        elif label == 'Image':
            densWindow.showImage = (not densWindow.showImage)
        elif label == 'Raw Bkg':
            densWindow.use_raw = (not densWindow.use_raw)

        densWindow.redraw ()




############################
# Task main()
############################

if __name__ == '__main__':

    # Process commandline args.
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hdvl:c:r:d:s:S:B:F:C:R:Z:",
            ["help","debug","verbose","load","catalog=", "ra=","dec=","sz=",
             "small_k=", "big_k=", "floor=", "ceiling=", "reticle=", "zoomsz="])

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

       elif opt in ("-c", "--catalog"):		# Catalog to query
           catalog = arg
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

       elif opt in ("-R", "--reticle"):		# Set cursor size
           cursz  = float (arg)
       elif opt in ("-Z", "--zoomsz"):		# Set size of zoom box
           zoomsz = float (arg)


    # Get the data and compute the differential convolution density
    # histograms.
    if fname is not None:
        ra0, dec0, g0, gr0, data0 = readData ( fname )

    elif xpos is not None and ypos is not None:
        ra0, dec0, g0, gr0, data0 = getData ( catalog )

    ra_idx = ra0 < 361				# Initialize for all values
    dec_idx = dec0 < 91


    # Initialize the density map.  [FIXME:  need to turn this into object]
    dmap = densMap (ra0, dec0, small_k, big_k)

    raw, extent, delta, clipped = dwarf_filter (dmap.ra, dmap.dec,
            ra_idx, dec_idx, fwhm_small=dmap.small_k, fwhm_big=dmap.big_k)
    dmap.setData (ra0, dec0)
    dmap.setHist (raw, extent, delta, clipped)


    # Create the subplot windows.
    fig = plt.figure ('Snow White...For finding small companions unexpectedly',
        figsize=(12,8))

    axdens = fig.add_subplot (121, aspect='equal')
    axzoom = fig.add_subplot (222, xlim=(xpos-xsize,xpos+xsize), 
        ylim=(ypos-ysize,ypos+ysize))
    axcmd  = plt.axes ([0.55, 0.1, 0.2, 0.35])
    axhist = plt.axes ([0.77, 0.1, 0.2, 0.35])

    # Create the density convolution subplot window.
    densWindow = densPlot (axdens)
    densWindow.setLabels ('RA (deg)', 'Dec (deg)')
    densWindow.setData (ra0, dec0)
    densWindow.redraw ()

    # Create the zoomed-region subplot window.
    zoomWindow = zoomPlot (axzoom)
    zoomWindow.setLabels ("RA", "Dec")
    zoomWindow.setData (ra0, dec0)
    zoomWindow.draw ()
    zoomWindow.drawContour (clipped, extent)

    # Create a plot of the CMD centered on the initial position.
    cmdWindow = cmdPlot (axcmd)
    cmdWindow.setLabels ("(g-r)", "g")
    cmdWindow.setData (gr0, g0)

    # Create the histogram subplot window.
    histWindow = histPlot (axhist, delta)
 
    # Create some task controls, explicitly positioned over the display.
    controls = controlFrame ()

    # Finalize the plot and display
    plt.tight_layout (pad=2.0)
    plt.show ()

