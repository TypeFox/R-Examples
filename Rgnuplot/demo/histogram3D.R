
# Initialize the gnuplot handle
h1 <- Gpinit()

# set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)

# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)

# filename for a temporary path+file
tmpfile <- tempfile()



# load the histogram3D script, set the filenames and number of rows and columns
Gpcmd(h1, "reset\nfile = \"immigr5.dat\"\ncylinder = \"" %s% tmpfile %s% "\"\ncol = 4\nrow = 5\nload \"histogram3D.gnu\"")

# pause R and gnuplot
Gppause()

# close gnuplot handle
h1 <- Gpclose(h1) 
