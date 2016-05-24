
# saving a plot to a PNG file (from gnuplot_i examples) Initialize the gnuplot handle
h1 <- Gpinit()
# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1, "set terminal png")
Gpcmd(h1, "set output \"sine.png\"")
GpplotEquation(h1, "sin(x)", "Sine wave")
# close gnuplot handles
g <- Gpclose(h1)
 
