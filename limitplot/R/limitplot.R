#
# Limitplot functions -- main constructor
#

limitplot <- function(..., lod, CI = 95, ratio = 1/25, shape = 1, size = 1, col = "black", main = "", xlab = "", ylab = "", names = "", axis = 5, stack = 5, jitterwidth = 0.2, jittershape = 1, jittersize = 1, jittercol = "black", log = "", blod = 1/2) {
  # Check LOD value
  if(lod<=0) {
    stop("Limit of detection (LOD) must be above 0.\n");
  }

  # Check values
  if(min(c(...))<=0) {
    stop("Variable values must be above 0.\n");
  }
 
  # Variable count
  vnum <- length(list(...));  

  # Set up the arrays with the variable data
  xa <- rep(seq(1:vnum), times = as.numeric(summary(list(...))[1:vnum]));
  ya <- c(...);

  if(log == "y") {
    ya <- log(ya);
    lod <- log(lod);
  }

  # Data frame with properly formatted data
  pl <- data.frame(xi = xa ,yi = ya);

  # Plot
  plotInit(pl, lod, vnum, stack, ratio, axis, main, xlab, ylab, jitterwidth, jittershape, jittersize, jittercol, names, log);
  plotBox(pl, lod, vnum, CI, blod);
  plotStack(pl, lod, vnum, stack, ratio, shape, size, col);
}
