
		CMD Magnifier -- Tools and Utilities

This repo contains the following files:

    cmdmag.py		-- Main interactive tool (i.e. Snow White)
    fieldplot.py	-- batch plotting application
    pal5.mov		-- demo movie of CMDMAG in action


      CMDMAG is a toy program to demo a query of the Data Lab TAP service 
and visualize the object density and associated CMD.  The user moves the 
cursor in the left-hand plot, the smaller plots track to show the catalog 
points zoomed in around the cursor, and the g-r vs g CMD for the points 
around the cursor, and a histogram of the differential values.

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

       -h, --help               # Help
       -d, --debug              # Debug
       -v, --verbose            # Verbose

       -c, --catalog            # Catalog to query
       -l, --load               # Load the FITS bintable
       -p, --ra                 # Set RA position
       -r, --dec                # Set DEC position
       -s, --sz                 # Set field size

       -S, --small_k            # Set small kernel size
       -B, --big_k              # Set big kernel size
       -F, --floor              # Set clipping floor (sigma)
       -C, --ceiling            # Set clipping ceiling (sigma)

       -R, --reticle            # Set cursor size
       -Z, --zoomsz             # Set size of zoom box

  Density-plot Commands:
        +                       # zoom in on field
        -                       # zoom out on field
        h / left-arrow          # pan left (get new data)
        j / down-arrow          # pan down (get new data)
        k / up-arrow            # pan up (get new data)
        l / right-arrow         # pan right (get new data)

        Shift                   # freeze motion tracking in density plot
        Shift-release           # enable motion tracking in density plot

  Zoom-plot Commands:
        < / >                   # change sample cursor by a factor of 2
        , / .                   # change sample cursor by a factor of 1.5
        +                       # increase zoom factor
        -                       # decrease zoom factor


Additional work to pre-fetch/cache data, or use of different density tools
will help speed and precision but is TBD.
