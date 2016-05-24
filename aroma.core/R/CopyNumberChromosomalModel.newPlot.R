setMethodS3("newPlot", "CopyNumberChromosomalModel", function(this, xlim, ylim=c(-1,1)*2.5, xlab="Physical position", ylab="Relative copy number", flavor=c("ce", "minimal"), xmargin=c(50,50), resScale=1, unit=6, yaxt="s", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'flavor':
  flavor <- match.arg(flavor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Reset graphical parameters when done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  xScale <- 1/(10^unit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plotting flavor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  args <- list();
  if (flavor == "glad") {
    par(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3));
  } else if (flavor == "ce") {
    # Margins in pixels-to-inches
    args <- c(args, list(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3)));

    # x-axis
    args <- c(args, list(xaxs="i"));

    # Set the horizontal margins to 'xmargin'.
    dim <- getDeviceResolution(resScale) * par("din");
    plt <- par("plt");    
    plt[1:2] <- c(xmargin[1], dim[1]-xmargin[2]) / dim[1];

    args <- c(args, "plt"=plt);
  } else if (flavor == "minimal") {
    # No axis
    yaxt <- "n";

    # No margins
    args <- c(args, list(mar=c(2,0,0.5,0), mgp=c(2,0.6,0.3)));

    # x-axis
    args <- c(args, list(xaxs="i"));
  }

  opar <- par(args);

  # Create empty plot figure
  plot(NA, xlim=xlim, ylim=ylim, xaxt="n", yaxt=yaxt, xlab=xlab, ylab=ylab, bty="n");

  invisible(opar);
}, protected=TRUE) # newPlot()


############################################################################
# HISTORY:
# 2007-10-17
# o Renamed to CopyNumberChromosomalModel.
# 2007-09-04
# o Created specifically for CopyNumberSegmentationModel.
# o Updated.  Now independent of 'fit' object.
# 2007-08-23
# o Created from profileCGH.plotProfile2.R.
############################################################################
