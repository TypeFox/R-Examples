setMethodS3("writeRegions" ,"profileCGH", function(this, filename, path=NULL, append=FALSE, ...) {
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Argument 'append':
  append <- Arguments$getLogical(append);


  pv <- this$profileValues;

  hasUnits <- (!is.null(pv$chipType) && !is.null(pv$units));
  if (hasUnits) {
    chipType <- as.character(pv$chipType);
    chipType <- gsub("[,-]monocell$", "", chipType);
    unitNames <- character(nrow(pv));
    for (cc in unique(chipType)) {
      cdfFile <- AffymetrixCdfFile$findByChipType(cc);
      if (is.null(cdfFile))
        throw("Cannot located CDF file for chip type: ", cc);
      idxs <- which(chipType == cc);
      unitNames[idxs] <- .readCdfUnitNames(cdfFile, units=pv$units[idxs]);
    }
  }
  snpNames <- NULL;
  rsIds <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write header?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!append) {
    cat(file=pathname, "# Breakpoints from GLAD fit\n");
    header <- "Chromosome\tstart\tstop\tlength\tGNL\tSmoothing\tnbrOfSNPs";
    if (hasUnits) {
      header <- paste(header, "firstSNP\tlastSNP", sep="\t");
      if (!is.null(rsIds))
        header <- paste(header, "firstRsId\tlastRsId", sep="\t");
    }
    header <- paste(header, "\n", sep="");
    cat(file=pathname, header, append=TRUE);
  }

  # Extract data to be exported
  regions <- pv$Region;
  chromosome <- pv$Chromosome;
  posBase <- pv$PosBase;
  zoneGNL <- pv$ZoneGNL;
  smoothing <- pv$Smoothing;

  # Identify unique regions
  uRegions <- unique(regions);

  # Format string
  fmtstr <- "%s\t%.f\t%.f\t%.f\t%s\t%+.2f\t%d";
  if (hasUnits) {
    fmtstr <- paste(fmtstr, "%s\t%s", sep="\t");
    if (!is.null(rsIds))
      fmtstr <- paste(fmtstr, "%s\t%s", sep="\t");
  }
  fmtstr <- paste(fmtstr, "\n", sep="");

  # For each region
  for (rr in uRegions) {
    # Get the first and last position of each region
    idx <- which(rr == regions);
    idx <- idx[c(1,length(idx))];

    if (hasUnits) {
      # Get the SNP names
      snpNames <- unitNames[idx];

    # Get the rsIDs
#      rsIds <- allRsIds[subset[idx]];
    }

    # Chromosome
    cc <- chromosome[idx[1]];

    # Is it a loss, normal or a gain?
    gnl <- c("loss", "normal", "gain")[zoneGNL[idx[1]]+2];

    # Create one line record
    args <- list(fmtstr, cc,
                 posBase[idx[1]], posBase[idx[2]], diff(posBase[idx]),
                 gnl, smoothing[idx[1]], as.integer(diff(idx)+1));
    if (hasUnits) {
      args <- c(args, snpNames[1], snpNames[2]);
      if (!is.null(rsIds))
        args <- c(args, rsIds[1], rsIds[2]);
    }
    record <- do.call(sprintf, args);
    cat(record);

    # Append to file
    cat(file=pathname, record, append=TRUE);
  }

  invisible(pathname);
})
