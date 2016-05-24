setMethodS3("pointsRawCNs", "default", function(fit, pchSymbol=c(20, 4), unit=6, col=NULL, ...) {
  rawCns <- extractRawCopyNumbers(fit);

  if (is.null(col)) {
    # TO DO: Bring in colors into a generic framework. /HB 2007-09-04
    cols <- RColorBrewer::brewer.pal(5, "Dark2");
    col <- cols[5];
  }

  points(rawCns, xScale=1/10^unit, pch=pchSymbol[1], col=col);
}, protected=FALSE);



# Patch for plotProfile() of class profileCGH so that 'ylim' argument works.
# Added also par(cex=0.8) - see code.
setMethodS3("pointsRawCNs", "profileCGH", function(fit, variable="LogRatio", chromosome=NULL, GNL="ZoneGNL", colDAGLAD=NULL, pchSymbol=c(20, 4), unit=6, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':
  if (!"PosBase" %in% names(fit$profileValues))
    throw("Argument 'fit' does not contain a 'PosBase' field.");

  # Argument 'variable':
  if (!variable %in% names(fit$profileValues))
    throw("Argument 'variable' does not specify a known field: ", variable);

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- unique(fit$profileValues$Chromosome);
    if (length(chromosome) > 1) {
      throw("Argument 'chromosome' must not be NULL if 'fit' contains more than one chromosome: ", paste(chromosome, collapse=", "));
    }
  }
  if (length(chromosome) > 1) {
    throw("Argument 'chromosome' must not contain more than one chromosome: ", paste(chromosome, collapse=", "));
  }

  # Argument 'colDAGLAD':
  if (is.null(colDAGLAD)) {
    colDAGLAD <- RColorBrewer::brewer.pal(5, "Dark2");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data for the chromosome to be plotted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pv <- fit$profileValues;

  # Keep only data for the chromosome of interest
  keep <- (pv$Chromosome == chromosome);
  pv <- pv[keep,];

  # Keep only finite values based on the variable of interest
  keep <- is.finite(pv[[variable]]);
  pv <- pv[keep,];

  # Make sure the order of the values are increasing (not really needed)
  o <- order(pv$PosOrder);
  pv <- pv[o,];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Annotate gains, normals, and losses, as well as outliers?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (GNL %in% names(pv)) {
    # Setup color vector
    col <- rep(colDAGLAD[5], length(pv$PosOrder));
    gnl <- pv[GNL];
    col[gnl ==  -1] <- colDAGLAD[4];
    col[gnl ==   1] <- colDAGLAD[3];
    col[gnl ==   2] <- colDAGLAD[2];
    col[gnl == -10] <- colDAGLAD[1];

    # Setup pch vector
    pch <- rep(pchSymbol[1], length(pv$PosOrder));
    pch[pv$OutliersTot != 0] <- pchSymbol[2];
  } else {
    col <- par("col");
    pch <- 20;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot the data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Extract the data to plot
  y <- pv[, variable];
  xScale <- 1/(10^unit);
  x <- xScale*pv$PosBase;

  points(x=x, y=y, pch=pch, col=col, ...);
}) # pointsRawCNs()



# Patch for plotProfile() of class profileCGH so that 'ylim' argument works.
# Added also par(cex=0.8) - see code.
setMethodS3("plotRawCNs", "profileCGH", function(fit, chromosome=NULL, unit=0, xlim=NULL, ylim=c(-1,1)*2.5, xlab="Physical position", ylab="Relative copy-number", flavor=c("glad", "ce", "minimal"), xmargin=c(50,50), resScale=1, ..., add=FALSE) {
  requireWithMemory("GLAD") || throw("Package not loaded: GLAD"); # data("cytoband")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fit':
  if (!"PosBase" %in% names(fit$profileValues))
    throw("Argument 'fit' does not contain a 'PosBase' field.");

  # Argument 'chromosome':
  if (is.null(chromosome)) {
    chromosome <- unique(fit$profileValues$Chromosome);
    if (length(chromosome) > 1) {
      throw("Argument 'chromosome' must not be NULL if 'fit' contains more than one chromosome: ", paste(chromosome, collapse=", "));
    }
  }
  if (length(chromosome) > 1) {
    throw("Argument 'chromosome' must not contain more than one chromosome: ", paste(chromosome, collapse=", "));
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get chromosome lengths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load data
  # To please R CMD check on R v2.6.0
  cytoband <- NULL; rm(list="cytoband");
  data("cytoband", envir=sys.frame(sys.nframe()));  # Package 'GLAD'
  genomeInfo <- aggregate(cytoband$End, list(Chromosome=cytoband$Chromosome,
                          ChrNumeric=cytoband$ChrNumeric), max, na.rm=TRUE);
  names(genomeInfo) <- c("Chromosome", "ChrNumeric", "Length");
  genomeInfo$Chromosome <- as.character(genomeInfo$Chromosome);
  genomeInfo$ChrNumeric <- as.integer(as.character(genomeInfo$ChrNumeric));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plotting flavor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reset graphical parameters when done
  opar <- par(no.readonly=TRUE);
  on.exit(opar);

  xScale <- 1/(10^unit);

  if (flavor == "glad") {
    par(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3));
    axes <- TRUE;
    if (is.null(xlim))
      xlim <- c(0, xScale*genomeInfo$Length[chromosome]);
  } else if (flavor == "ce") {
    # Margins in pixels-to-inches

    par(mar=c(3,3,5,3)+0.1, mgp=c(2,0.6,0.3), xaxs="i");

    # Set the horizontal margins to 'xmargin'.
    dim <- getDeviceResolution(resScale) * par("din");
    plt <- par("plt");
    plt[1:2] <- c(xmargin[1], dim[1]-xmargin[2]) / dim[1];
    par("plt"=plt);

    axes <- TRUE;
    if (is.null(xlim))
      xlim <- c(0, xScale*genomeInfo$Length[chromosome]);
  } else if (flavor == "minimal") {
    # No margins
    	par(mar=c(2,0,0.5,0), mgp=c(2,0.6,0.3), xaxs="i");
    # No axis
    axes <- FALSE;
    # No cytobands
    plotband <- FALSE;
    # x-range
    if (is.null(xlim))
      xlim <- c(0, xScale*genomeInfo$Length[chromosome]);
  }

  # Create plot
  if (!add) {
    plot(NA, xaxt="n", xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=axes);
  }

  # Plot raw CN data
  pointsRawCNs(fit, chromosome=chromosome, unit=unit, ...);
}) # plotRawCNs()


############################################################################
# HISTORY:
# 2012-11-29
# o Dropped 'verbose' usage in pointsRawCNs() for default.
# 2007-09-04
# o Added a default pointsRawCNs().
# 2007-08-22
# o Created from plotProfile2().
############################################################################
