setMethodS3("lapplyInChunks", "numeric", function(idxs, fcn, ..., chunkSize=1, useNames=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'fcn':
  if (!is.function(fcn)) {
    throw("Argument 'fcn' is not a function: ", mode(fcn));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "lapplyInChunks()");

  n <- length(idxs);

  # Nothing todo?
  if (n == 0)
    return(list());

  remaining <- 1:n;
  chunk <- 1:chunkSize;
  nbrOfChunks <- ceiling(n / chunkSize);

  # Allocate return structure
  res <- vector("list", n);
  names <- NULL;

  verbose && cat(verbose, "Number of elements per chunk: ", chunkSize);

  count <- 1;
  while (length(remaining) > 0) {
    verbose && enter(verbose, "Chunk #", count, " of ", nbrOfChunks);
    # Last chunk?
    if (length(remaining) < chunkSize)
      chunk <- 1:length(remaining);

    # Indices for this chunk
    ii <- remaining[chunk];

    verbose && cat(verbose, "Elements: ");
    verbose && str(verbose, idxs[ii]);

    # Apply the function
    resChunk <- fcn(idxs[ii], ...);
    res[ii] <- resChunk;

    if (useNames && is.character(names(resChunk))) {
      if (is.null(names))
        names <- vector("character", n);
      names[ii] <- names(resChunk);
    }

    # Not needed anymore
    resChunk <- ii <- NULL;

    # Next chunk
    remaining <- remaining[-chunk];
    count <- count + 1;

    verbose && exit(verbose);
  } # while(length(idxs) > 0)

  # Not needed anymore
  chunk <- NULL;

  # Add names?
  if (!is.null(names))
    names(res) <- names;

  verbose && exit(verbose);

  res;
}) # lapplyInChunks()


setMethodS3("splitInChunks", "numeric", function(idxs, chunkSize=1, ...) {
  n <- length(idxs);
  nbrOfChunks <- ceiling(n / chunkSize);
  res <- vector("list", nbrOfChunks);
  head <- seq_len(chunkSize);
  pos <- 1;
  while (length(idxs) > 0) {
    if (length(idxs) < chunkSize)
      head <- seq_along(idxs);
    res[[pos]] <- idxs[head];
    idxs <- idxs[-head];
    pos <- pos + 1;
  }
  res;
}) # splitInChunks()


############################################################################
# HISTORY:
# 2011-03-03
# o BUG FIX: lapplyInChunks(idxs) for numeric did not correctly handle
#   the case when length(idxs) == 0, because of a typo.
# 2007-08-17
# o Added splitInChunks().
# 2007-02-12
# o Created.
############################################################################
