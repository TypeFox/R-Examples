
#'Hello World!' - text on caption
# Initialize the gnuplot handle
h1 <- Gpinit()
# naming the axis
GpsetXlabel(h1, "x")
GpsetYlabel(h1, "y")
# set plot style to 'lines'
Gpsetstyle(h1, "lines")
# add a caption
Gpcmd(h1, "set title \"Hello World!\"")
# plot sin(x)
GpplotEquation(h1, "sin(x)", "")
# pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
