
# 'Hello World!' - using GpplotEquation
h1 <- Gpinit()  #Initialize the gnuplot handle
GpplotEquation(h1, "sin( x )", "Hello World!")
Gppause()  # pause R and gnuplot
h1 <- Gpclose(h1)  # close gnuplot handles 
