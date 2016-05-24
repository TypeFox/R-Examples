
# example of using GpLoadDemo and GpReadURL2string Initialize the gnuplot handle
h1 <- Gpinit()
# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
# load the file 'simple.dem' Gpcmd(h1, 'set terminal postscript eps color;set output 'simple.eps'\n' %s% GpURL2string('http://gnuplot.sourceforge.net/demo_svg/simple.1.gnu') %s%
# '\nset terminal X11;set output')
if (!file.exists("/usr/share/doc/gnuplot-doc/examples/simple.dem")) stop("Please install gnuplot-doc")
GpLoadDemo(h1, "/usr/share/doc/gnuplot-doc/examples/simple.dem")
# pause R and gnuplot
Gppause()
# example of GpReadURL2string Kuen's Surface
gpcode <- GpURL2string("http://gnuplot.sourceforge.net/demo/transparent_solids.2.gnu")
# send gnuplot script
Gpcmd(h1, gpcode)
# Gpcmd(h1, 'set terminal postscript eps color;set output 'KuensSurface.eps'\n' %s% gpcode) pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
