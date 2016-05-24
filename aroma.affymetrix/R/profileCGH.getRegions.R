setMethodS3("getRegions", "profileCGH", function(this, nbrOfSnps=c(1,Inf), smoothing=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'nbrOfSnps':
  if (length(nbrOfSnps) == 1)
    nbrOfSnps <- c(nbrOfSnps,Inf);

  # Argument 'smoothing':
  if (!is.null(smoothing) & !is.matrix(smoothing)) {
    smoothing <- matrix(smoothing, ncol=2, byrow=TRUE);
  }

  pv <- this$profileValues;
  stdvs <- this$SigmaC$Value;

  hasUnits <- (!is.null(pv$chipType) && !is.null(pv$units));
  if (hasUnits) {
    chipType <- as.character(pv$chipType);
    chipType <- gsub("[,-]monocell$", "", chipType);

    rsIds <- character(nrow(pv));
    unitNames <- character(nrow(pv));
    for (cc in unique(chipType)) {
      cdf <- AffymetrixCdfFile$byChipType(cc);

      idxs <- which(chipType == cc);
      unitNames[idxs] <- getUnitNames(cdf, units=pv$units[idxs]);
    }
  }

  rsIds <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate result table
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify unique regions
  uRegions <- unique(pv$Region);
  nbrOfRegions <- length(uRegions);

  # Columns
  colClasses <- c(Chromosome="character", start="integer", 
                   stop="integer", length="integer", nbrOfSnps="integer", 
                                       Smoothing="double", SNRtoZero="double");
  if (hasUnits) {
    colClasses <- c(colClasses, firstSnp="character", lastSnp="character");
    if (!is.null(rsIds))
      colClasses <- c(colClasses, firstRsId="character", lastRsId="character");
  }

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

    # Chromosome
    df[rr,"Chromosome"] <- pv$Chromosome[idx[1]];

    # (start, stop, length)
    df[rr,c("start", "stop")] <- as.integer(pv$PosBase[idx]);
    df[rr,"length"] <- as.integer(diff(pv$PosBase[idx]));

    # Number of SNPs
    df[rr,"nbrOfSnps"] <- as.integer(diff(idx)+1);

    # Smoothing
    df[rr,"Smoothing"] <- pv$Smoothing[idx[1]];

    # Signal-to-noise ratio
    df[rr,"SNRtoZero"] <- abs(df[rr,"Smoothing"]) / stdvs;

    # Gain, normal, or loss?
#    levels <- c("loss", "normal", "gain");
#    levels <- c("-", "0", "+");
#    df[rr,"GNL"] <- levels[pv$ZoneGNL[idx[1]]+2];

    if (hasUnits) {
      # Get the SNP names
      df[rr,c("firstSnp", "lastSnp")] <- unitNames[idx];

      # Get the rsIds 
#      if (!is.null(rsIds))
#        df[rr,c("firstRsId", "lastRsId")] <- allRsIds[subset[idx]];
# 'allRsIds', 'subset'?!? /HB 2007-06-11
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filter
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Gain, normal, or loss
#  if (!is.null(gnl)) {
#    keep <- (df$GNL %in% gnl);
#    df <- df[keep,];
#  }

  # Number of SNPs
  if (!is.null(nbrOfSnps)) {
    keep <- (nbrOfSnps[1] <= df$nbrOfSnps & df$nbrOfSnps <= nbrOfSnps[2]);
    df <- df[keep,];
  }

  # Smoothing regions
  if (!is.null(smoothing)) {
    keep <- rep(FALSE, nrow(df));
    for (kk in seq_len(nrow(smoothing))) {
      range <- smoothing[kk,];
      keep <- keep | (range[1] <= df$Smoothing & df$Smoothing <= range[2]);
    }
    df <- df[keep,];
  }

  df;
}, private=TRUE)

##############################################################################
# HISTORY:
# 2007-06-11
# o RS IDs are currently not returned/supported.
# 2006-12-11
# o Added SNR column.
# 2006-11-22
# o Created from writeRegions.profileCGH.R from 2006-10-31.
##############################################################################
