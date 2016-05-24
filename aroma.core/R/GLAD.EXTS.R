setMethodS3("extractCopyNumberRegions", "profileCGH", function(object, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pv <- object$profileValues;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate result table
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify unique regions
  uRegions <- unique(pv$Region);
  nbrOfRegions <- length(uRegions);

  # Columns
  colClasses <- c(chromosome="character", start="integer",
                  stop="integer", mean="double", nbrOfLoci="integer",
                  call="character");
  df <- dataFrame(colClasses, nrow=nbrOfRegions);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract each region
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (rr in seq_along(uRegions)) {
    # Get the region ID
    region <- uRegions[rr];

    # Get the first and last position of each region
    idx <- which(region == pv$Region);
    idx <- idx[c(1,length(idx))];
    idx1 <- idx[1];

    # Chromosome
    df[rr,"chromosome"] <- pv$Chromosome[idx1];

    # (start, stop, length)
    df[rr,c("start", "stop")] <- as.integer(pv$PosBase[idx]);

    # Number of SNPs
    df[rr,"nbrOfLoci"] <- as.integer(diff(idx)+1);

    # Smoothing
    df[rr,"mean"] <- pv$Smoothing[idx1];

    # Call
    df[rr,"call"] <- c("loss", "neutral", "gain")[pv$ZoneGNL[idx1]+2];
  }

  CopyNumberRegions(
    chromosome=df$chromosome,
    start=df$start,
    stop=df$stop,
    mean=df$mean,
    count=df$nbrOfLoci,
    call=df$call
  );
}) # extractCopyNumberRegions()


setMethodS3("extractRawCopyNumbers", "profileCGH", function(object, ...) {
  pv <- object$profileValues;
  chromosome <- unique(pv$Chromosome);
  chromosome <- Arguments$getIndex(chromosome);
  RawCopyNumbers(cn=pv$LogRatio, x=pv$PosBase, chromosome=chromosome);
})


setMethodS3("drawCnRegions", "profileCGH", function(this, ...) {
  cnr <- extractCopyNumberRegions(this, ...);
  drawLevels(cnr, ...);
})


# Patch for plotProfile() of class profileCGH so that 'ylim' argument works.
# Added also par(cex=0.8) - see code.
setMethodS3("drawCytoband", "profileCGH", function(fit, chromosome=NULL, cytobandLabels=TRUE, colCytoBand=c("white", "darkblue"), colCentro="red", unit=6, ...) {
  requireWithMemory("GLAD") || throw("Package not loaded: GLAD");  # data("cytoband")

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


  xScale <- 1/(10^unit);


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
  # Get the cytoband details for the chromosome of interest
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Drop column 'Chromosome'
  ## Gives a NOTE in R CMD check R v2.6.0, which is nothing, but we'll
  ## use a workaround to get a clean result. /HB 2007-06-12
  Chromosome <- NULL; rm(list="Chromosome"); # dummy
  cytobandNew <- subset(cytoband, select=-Chromosome);
  cytobandNew <- merge(LabelChr, cytobandNew, by.x="Chromosome",
                                                        by.y="ChrNumeric");
  # Rescale x positions according to units
  cytobandNew$Start <- xScale*cytobandNew$Start;
  cytobandNew$End <- xScale*cytobandNew$End;

  # Where should the cytoband be added and how wide should it be?
  usr <- par("usr");
  dy <- diff(usr[3:4]);

  drawCytoband2(cytobandNew, chromosome=chromosome,
    labels=cytobandLabels, y=usr[4]+0.02*dy, height=0.03*dy,
    colCytoBand=colCytoBand, colCentro=colCentro);
}, private=TRUE) # drawCytoband()




############################################################################
# HISTORY:
# 2010-02-19
# o Moved drawCytoband2() to its own file, because it no longer requires
#   the GLAD package.
# 2009-05-14
# o Moved extractRawCopyNumbers() for profileCGH from aroma.affymetrix.
# o Moved extractCopyNumberRegions() for profileCGH from aroma.affymetrix.
# 2009-05-10
# o Moved to aroma.core v1.0.6.  Source files: profileCGH.drawCnRegions.R
#   and profileCGH.drawCytoband.R.
# 2008-05-21
# o Now extractRawCopyNumbers() adds 'chromosome' to the returned object.
# 2007-09-04
# o Now data("cytoband") is loaded to the local environment.
# 2007-08-22
# o Update plotProfile2() to utilizes drawCnRegions().
# 2007-08-20
# o Added drawCnRegions().
# 2007-06-11
# o Added explicit call to GLAD::myPalette() to please R CMD check R v2.6.0.
# 2007-01-03
# o Made the highlighting "arrow" for the centromere smaller.
# 2006-12-20
# o It is now possible to specify 'xlim' as well as 'ylim'.
# o Reimplemented, because the cytoband was not displayed correctly.
############################################################################
