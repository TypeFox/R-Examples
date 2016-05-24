
#'Hello World!' - loading a script directly from gnuplot
# Initialize the gnuplot handle
h1 <- Gpinit()
# set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)
# load the script
Gpcmd(h1, "load \"helloworld.gnu\"")
# pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
