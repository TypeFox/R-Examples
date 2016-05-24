# Patch for plotProfile() of class profileCGH so that 'ylim' argument works.
# Added also par(cex=0.8) - see code.
setMethodS3("plotProfile2", "profileCGH", function(fit, variable="LogRatio", chromosome=NULL, Smoothing="Smoothing", GNL="ZoneGNL", Bkp=FALSE, cytobandLabels=TRUE, plotband=TRUE, unit=0, colDAGLAD=NULL, pchSymbol=c(20, 4), colCytoBand=c("white", "darkblue"), colCentro="red", xlim=NULL, ylim=c(-1,1)*2.5, xlab="Physical position", ylab=variable, flavor=c("glad", "ce", "minimal"), xmargin=c(50,50), resScale=1, ...) {
  requireWithMemory("GLAD") || throw("Package not loaded: GLAD"); # data("cytoband")

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

  # Argument 'Smoothing':
  if (!is.null(Smoothing)) {
    if (!Smoothing %in% names(fit$profileValues)) {
      cat("Warning in plotProfile.profileCGH:", Smoothing, " is not available");
    }
  }

  # Argument 'Bkp':
  if (Bkp) {
    if (!"Breakpoints" %in% names(fit$profileValues))
      throw("Cannot plot breakpoints: No data available.")
  }

  # Argument 'colDAGLAD':
  if (is.null(colDAGLAD)) {
    colDAGLAD <- RColorBrewer::brewer.pal(5, "Dark2");
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reset graphical parameters when done
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  xScale <- 1/(10^unit);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Keep only data to be plotted
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pv <- fit$profileValues;

  # Keep only data for the chromosome of interest
  keep <- (pv$Chromosome == chromosome);
  pv <- pv[keep,];

  # Keep only finite values based on the variable of interest
  keep <- is.finite(pv[[variable]]);
  pv <- pv[keep,];

  # Convert the chromosome names to chromosome indices
  pv$Chromosome <- GLAD::ChrNumeric(pv$Chromosome);

  # Make sure the order of the values are increasing
  o <- order(pv$PosOrder);
  pv <- pv[o,];


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get chromosome lengths
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load data
  # To please R CMD check on R v2.6.0
  cytoband <- NULL; rm(list="cytoband");
  data("cytoband", envir=sys.frame(sys.nframe()));  # Package 'GLAD'
  genomeInfo <- aggregate(cytoband$End,
    by=list(Chromosome=cytoband$Chromosome, ChrNumeric=cytoband$ChrNumeric),
    FUN=max, na.rm=TRUE);
  names(genomeInfo) <- c("Chromosome", "ChrNumeric", "Length");
  genomeInfo$Chromosome <- as.character(genomeInfo$Chromosome);
  genomeInfo$ChrNumeric <- as.integer(as.character(genomeInfo$ChrNumeric));

  LabelChr <- data.frame(Chromosome=chromosome);
  LabelChr <- merge(LabelChr, genomeInfo[, c("ChrNumeric", "Length")],
                         by.x="Chromosome", by.y="ChrNumeric", all.x=TRUE);

  LabelChr$Length <- 0;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the plot data with cytoband information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pv <- merge(pv, LabelChr, by="Chromosome");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plotting flavor
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (flavor == "glad") {
  } else if (flavor == "ce") {
  } else if (flavor == "minimal") {
    # No cytobands
    plotband <- FALSE;
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create empty plot figure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  newPlot(fit, unit=unit, xlim=xlim, ylim=ylim, flavor=flavor, ...);
#  plot(NA, xlim=xlim, ylim=ylim, xaxt="n", xlab=xlab, ylab=ylab, bty="n");


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

    colSmoothing <- "black";
  } else {
    col <- par("col");
    pch <- 20;
    colSmoothing <- "red";
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot main data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the data to plot
  y <- pv[, variable];
  x <- xScale*pv$PosBase;

#   plot(x=x, y=y, pch=pch, col=col, xlim=xlim, ylim=ylim, xaxt="n", xlab=xlab, ylab=ylab, bty="n");

#  points(x=x, y=y, pch=pch, col=col, ...);
  pointsRawCNs(fit, unit=unit, ...);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot cytobands?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (plotband) {
    drawCytoband(fit, chromosome=chromosome, cytobandLabels=TRUE,
           colCytoBand=colCytoBand, colCentro=colCentro, unit=unit);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot break points?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (Bkp) {
    if (is.data.frame(fit$BkpInfo)) {
      fit$BkpInfo <- merge(fit$BkpInfo, LabelChr, by="Chromosome");
      fit$BkpInfo$NewPosBase <- fit$BkpInfo$PosBase + fit$BkpInfo$Length;
      fit$BkpInfo$NewPosBase <- xScale*fit$BkpInfo$NewPosBase;
      abline(v=fit$BkpInfo$NewPosBase + 0.5, col="red", lty=2);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot smoothing values?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(Smoothing)) {
    drawCnRegions(fit, ..., col=colSmoothing, xScale=xScale);
  }
}, private=TRUE) # plotProfile2()


############################################################################
# HISTORY:
# 2010-12-07
# o plotProfile2() for profileCGH now utilizing requireWithMemory()
#   to decrease the annoyances for users if GLAD fails to load.
# 2007-09-04
# o Now data("cytoband") is loaded to the local environment.
# 2007-08-22
# o Update plotProfile2() to utilizes drawCnRegions().
# 2007-06-11
# o Added explicit call to GLAD::myPalette() to please R CMD check R v2.6.0.
# 2007-01-03
# o Made the highlighting "arrow" for the centromere smaller.
# 2006-12-20
# o It is now possible to specify 'xlim' as well as 'ylim'.
# o Reimplemented, because the cytoband was not displayed correctly.
############################################################################
