# @author "MR"
setMethodS3("importFromBpMap", "AromaCellMatchScoreFile", function(this, srcPathname, rows=NULL, ..., verbose=TRUE) {
  # Argument 'srcPathname':
  srcPathname <- Arguments$getReadablePathname(srcPathname);

  # Argument 'rows':
  if(is.null(rows))
    stop("Must provide the chip dimensions: 'rows' argument is NULL");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Importing match scores from BPMAP file");

  verbose && enter(verbose, "Reading BPMAP file");
  verbose && cat(verbose, "Pathname: ", srcPathname);
  bps <- .readBpmap(srcPathname, readMatchScore=TRUE);
  verbose && exit(verbose);

  verbose && enter(verbose, "Saving scores to ACM file");
  verbose && cat(verbose, "Pathname: ", getPathname(this));

  for(kk in seq_len(length(bps))) {
    bp <- bps[[kk]];
    verbose && enter(verbose, "Updating ", bp$seqInfo$fullname[1]);

    ms <- round(bp$matchscore*1e6);
    # which() is faster than which()
    w <- which(ms >= 1 & ms <= 10);
    if(length(w) > 0) {
      cells <- bp$pmy[w]*rows + bp$pmx[w] + 1;
      updateMatchScores(this, cells=cells, scores=ms[w]);
    }
    # Not needed anymore
    ms <- w <- NULL;

    verbose && exit(verbose);
  } # for (kk ...)
  verbose && exit(verbose);

  invisible(this);
})


############################################################################
# HISTORY:
# 2008-10-xx [MR]
# o Created.
############################################################################
