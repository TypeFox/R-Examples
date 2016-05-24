setMethodS3("plotAllelePairs", "AllelicCrosstalkCalibration", function(this, array, pairs=NULL, what=c("input", "output"), ..., pch=".", cex=1, lty=3, xlim=c(-500,2^16), ylim=xlim, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'array':
  array <- Arguments$getIndex(array, range=c(1, Inf));

  # Argument 'pairs':
  if (is.null(pairs)) {
  } else {
    pairs <- Arguments$getIndices(pairs, max=999);
  }

  # Argument 'what':
  what <- match.arg(what);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get data set of interest
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (what == "input") {
    ds <- getInputDataSet(this);
  } else if (what == "output") {
    if (!isDone(this)) {
      throw("Cannot plot allele-pair signals for output data set. Data is still not calibrated.");
    }
    ds <- getOutputDataSet(this);
  }

  # Get the array of interest
  df <- ds[[array]];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the sets of cell-index pairs of interest
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setsOfProbes <- getSetsOfProbes(this, verbose=less(verbose, 5));
  setsOfProbes <- setsOfProbes$snps;
  nbrOfPairs <- length(setsOfProbes);
  if (is.null(pairs)) {
    pairs <- seq_len(nbrOfPairs);
  } else {
    pairs <- Arguments$getIndices(pairs, max=nbrOfPairs);
    nbrOfPairs <- length(pairs);
  }

  # Sanity check
  if (nbrOfPairs < 1) {
    throw("No allele pairs to plot");
  }

  # Extract cell pairs of interest
  setsOfProbes <- setsOfProbes[pairs];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set up the plot grid
  layout(matrix(seq_len(nbrOfPairs), nrow=floor(sqrt(nbrOfPairs))));
  par(mar=c(3,3,1,1)+0.1, mgp=c(1.8,0.5,0), oma=c(0,0,3,0));
  for (kk in seq_len(nbrOfPairs)) {
    cells <- setsOfProbes[[kk]];
    name <- names(setsOfProbes)[kk];
    y <- extractMatrix(df, cells=cells);
    dim(y) <- dim(cells);
    colnames(y) <- colnames(cells);
    xlab <- substitute(y[nn], list(nn=colnames(y)[1]));
    ylab <- substitute(y[nn], list(nn=colnames(y)[2]));
    plot(NA, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
    if (kk == 1)
      mtext(getFullName(df), outer=TRUE, line=1);
    stext(side=3, pos=1, sprintf("%s (n=%d)", name, nrow(y)));
    x1 <- par("usr")[2];
    y1 <- par("usr")[4];
    lines(x=c(0,x1), y=c(0,y1), lty=3, lwd=2, col="#999999");
    lines(x=c(0,0), y=c(0,y1), lty=3, lwd=2, col="#999999");
    lines(x=c(0,x1), y=c(0,0), lty=3, lwd=2, col="#999999");
    points(y, pch=pch, cex=cex, ...);
  } # for (kk ...)
  invisible(list(setsOfProbes=setsOfProbes, df=df));
}, protected=TRUE) # plotAllelePairs()



############################################################################
# HISTORY:
# 2008-12-10
# o Added plotAllelePairs() for AllelicCrosstalkCalibration.
# o Created.
############################################################################
