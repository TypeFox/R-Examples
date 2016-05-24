setMethodS3("extractLocusLevelC1C2", "PairedPSCBS", function(fit, ...) {
  data <- getLocusData(fit);
  C <- data$CT;
  rho <- data$rho;
  C1 <- 1/2*(1-rho)*C;
  C2 <- C-C1;
  data.frame(C1=C1, C2=C2);
}, private=TRUE) # extractLocusLevelC1C2()


setMethodS3("extractLocusLevelTCN", "PairedPSCBS", function(fit, ...) {
  data <- getLocusData(fit);
  C <- data$CT;
}, private=TRUE) # extractLocusLevelTCN()



setMethodS3("extractDhSegment", "PairedPSCBS", function(fit, idx, what=c("hets", "SNPs", "all"), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'what':
  what <- match.arg(what);


  segs <- getSegments(fit, splitters=TRUE);
  stopifnot(!is.null(segs)); 
  nbrOfSegments <- nrow(segs);

  # Argument 'idx':
  idx <- Arguments$getIndex(idx, max=nbrOfSegments);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 



  verbose && enter(verbose, "Extracting a specific DH segment");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Extract the data and segmentation results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  data <- getLocusData(fit);
  stopifnot(!is.null(data));

  segs <- getSegments(fit, splitters=TRUE);
  stopifnot(!is.null(segs));

  verbose && enter(verbose, "Subsetting segment");
  # Subset the region-level data
  seg <- segs[idx,,drop=FALSE];

  isDivider <- all(is.na(seg));
  if (isDivider) {
    verbose && cat("Cannot extract DH segment. Not a valid segment: ", idx);
    verbose && exit(verbose);
    return(NULL);
  }

  verbose && print(verbose, seg);
  verbose && cat(verbose, "Number of TCN markers: ", sum(seg[["tcnNbrOfLoci"]], na.rm=TRUE));
  verbose && exit(verbose);

  verbose && enter(verbose, "Subsetting data");
  units <- seq(length=nrow(data));

  # Keep only chromosome of interest
  chr <- as.numeric(seg[,"chromosome"]);
  if (!is.na(chr)) {
    keep <- which(data$chromosome == chr);
    units <- units[keep];
    data <- data[keep,];
  }

  # Keep only loci within the segment
  xRange <- as.numeric(seg[,c("dhStart", "dhEnd")]);
  keep <- which(xRange[1] <= data$x & data$x <= xRange[2]);
  units <- units[keep];
  data <- data[keep,];

  muN <- data$muN;
  isSnp <- is.finite(muN);

  # Keep only SNPs?
  if (is.element(what, c("SNPs", "hets"))) {
    keep <- which(isSnp);
    units <- units[keep];
    data <- data[keep,];
  }

  # Keep only heterozygous SNPs?
  if (what == "hets") {
    isHet <- (muN == 1/2);
    keep <- which(isHet);
    units <- units[keep];
    data <- data[keep,];
  }
  verbose && exit(verbose);

  n <- nrow(data);
  verbose && cat(verbose, "Number of loci in DH segment: ", n);

  # Special case?
  listOfDhLociNotPartOfSegment <- fit$listOfDhLociNotPartOfSegment;
  if (!is.null(listOfDhLociNotPartOfSegment)) {
    tcnId <- seg[,"tcnId"];
    dhId <- seg[,"dhId"];
    dhLociNotPartOfSegment <- listOfDhLociNotPartOfSegment[[tcnId]];
    if (!is.null(dhLociNotPartOfSegment)) {
      lociToExclude <- dhLociNotPartOfSegment[[dhId]];
      verbose && cat(verbose, "Excluding loci that belongs to a flanking segment: ", length(lociToExclude));
      drop <- match(lociToExclude, units);
      units <- units[-drop];
      data <- data[-drop,];
      n <- nrow(data);
    }
  }

  verbose && cat(verbose, "Number of units: ", n);
  verbose && cat(verbose, "Number of TCN markers: ", seg[,"tcnNbrOfLoci"]);

  # Sanity check
  if (what == "hets" && n > 0) stopifnot(n == seg[,"dhNbrOfLoci"]);

  fitS <- fit;
  fitS$data <- data;
  fitS$output <- seg;

  verbose && exit(verbose);

  fitS;
}, protected=TRUE) # extractDhSegment()


############################################################################
# HISTORY:
# 2012-02-24
# o Added extractDhSegment() for PairedPSCBS, which was copied "as is"
#   from the aroma.cn package.  The below history has been updated to
#   document changes in this method too.
# 2012-02-23
# o Made extractDhSegment() protected. 
# 2011-10-08
# o ROBUSTIFICATION: Uses drop=FALSE in mergeTwoSegments() for PairedPSCBS.
# 2010-10-26 [HB]
# o Added extractDhSegment() for PairedPSCBS.
# 2011-10-02
# o DOCUMENTATION: Added Rdoc help to mergeTwoSegments() & dropByRegions().
# o Added verbose statements to the above to functions.
# 2011-06-14
# o Updated code to recognize new column names.
# 2011-01-18
# o BUG FIX: Fields 'tcnSegRows' and 'dhSegRows' were not updated by
#   mergeTwoSegments() for PairedPSCBS.
# 2011-01-14
# o Moved extractByRegions() and estimateStdDevForHeterozygousBAF() to
#   psCBS v0.9.36.
# o Now extractByRegions() utilizes the 'segRows' field.
# o Added estimateStdDevForHeterozygousBAF().
# 2011-01-12
# o Added updateMeans() for PairedPSCBS.
# o Added dropByRegions().
# o Added extractByRegions() and extractByRegion().
# o Created.
############################################################################
