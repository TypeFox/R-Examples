# @author "MR"
 setMethodS3("importFromBpMap", "AromaCellSequenceFile", function(this, srcPathname, rows=NULL, ..., verbose=TRUE) {
  # Argument 'srcPathname':
  srcPathname <- Arguments$getReadablePathname(srcPathname);

  # Argument 'rows':
  # HB: Can't we infer this from the BPMAP file? ...or have an option
  # to look it up in matching CDF.
  if (is.null(rows)) {
    stop("Must provide the chip dimensions: 'rows' argument is NULL");
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing cell sequences from BPMAP file");

  verbose && enter(verbose, "Reading BPMAP file");
  verbose && cat(verbose, "Pathname: ", srcPathname);
  bps <- .readBpmap(srcPathname);
  verbose && exit(verbose);

  verbose && enter(verbose, "Saving to ACS file");
  verbose && cat(verbose, "Pathname: ", getPathname(this));
  for(kk in seq_len(length(bps))) {
    bp <- bps[[kk]];
    verbose && enter(verbose,"Updating ", bp$seqInfo$fullname[1]);

    cells <- bp$pmy*rows + bp$pmx + 1;
    seqs <- bp$probeseq;
    updateSequences(this, cells=cells, seqs=seqs);

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  verbose && exit(verbose);
  invisible(this);
})

############################################################################
# HISTORY:
# 2008-10-xx [MR]
# o Created.
############################################################################
