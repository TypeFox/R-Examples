local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

# **********************************************
# R code for the Lissajous example
# **********************************************

# Calculate and draw a Lissajous figure ("Plot" button)

drawLiss <- function() {
  getWinVal(scope="L");
  ti <- 2*pi*(0:k)/k;
  x <- sin(2*pi*m*ti);
  y <- sin(2*pi*(n*ti+phi));
  resetGraph();
  plot(x,y,type=ptype);
  invisible(NULL); }

# Load PBS Modelling and initialize the GUI

require(PBSmodelling); createWin("LissFigWin.txt")

}) # end local scope
