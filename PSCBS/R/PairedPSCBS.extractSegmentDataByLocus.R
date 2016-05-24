setMethodS3("extractSegmentDataByLocus", "PairedPSCBS", function(fit, fields=NULL, ..., verbose=FALSE) {
  # Extract data
  segs <- getSegments(fit, splitters=TRUE);
  data <- getLocusData(fit, ...);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fields':
  if (!is.null(fields)) {
    fields <- Arguments$getCharacters(fields);
    unknown <- fields[!is.element(fields, colnames(segs))];
    if (length(unknown) > 0L) {
      throw("Unknown segment fields: ", paste(sQuote(unknown), collapse=", "));
    }
  } else {
    fields <- colnames(segs);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting segment data by locus");

  # Extract segment fields
  chromosome <- data$chromosome;
  x <- data$x;
  y <- data[,3L];
  segs <- segs[,fields];
  nbrOfLoci <- nrow(data);

  verbose && printf(verbose, "Segment fields: [%d] %s\n", length(fields), paste(sQuote(fields), collapse=", "));
  verbose && cat(verbose, "Number of loci: ", nbrOfLoci);

  # Allocate segment fields at the locus level
  dataL <- matrix(NA, nrow=nbrOfLoci, ncol=length(fields));
  colnames(dataL) <- fields;
  dataL <- as.data.frame(dataL);

  verbose && cat(verbose, "Allocated results:");
  verbose && str(verbose, dataL);

  verbose && enter(verbose, "Extracting segment by segment");

  # For each segment...
  for (ss in seq(length=nrow(segs))) {
    verbose && enter(verbose, sprintf("Segment %d of %d", ss, nrow(segs)));
    seg <- segs[ss,];
    idxs <- which(chromosome == seg$chromosome &
                  seg$tcnStart <= x & x <= seg$tcnEnd);
    idxs <- Arguments$getIndices(idxs, max=nbrOfLoci);
    verbose && cat(verbose, "Number of loci in segment: ", length(idxs));
    # Sanity check
##    stopifnot(length(idxs) == seg$tcnNbrOfLoci);

    segsSS <- seg[fields];
    verbose && cat(verbose, "Segment data extracted:");
    verbose && print(verbose, segsSS);

    for (field in fields) {
      dataL[idxs,field] <- segsSS[[field]];
    }
    verbose && exit(verbose);
  } # for (ss ...)
  verbose && exit(verbose);

  # The calls for loci that have missing annotations or observations,
  # should also be missing, i.e. NA.
  nok <- (is.na(chromosome) | is.na(x) | is.na(y));
  dataL[nok,] <- NA;

  # Sanity check
  stopifnot(nrow(dataL) == nbrOfLoci);
  stopifnot(ncol(dataL) == length(fields));

  verbose && exit(verbose);

  dataL;
}, protected=TRUE) # extractSegmentDataByLocus()



##############################################################################
# HISTORY
# 2013-10-27
# o Added extractSegmentDataByLocus() for PairedPSCBS adopted from
#   extractCallsByLocus() of 'PSCBS'.
##############################################################################
