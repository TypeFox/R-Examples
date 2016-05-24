# @RdocClass "AromaCellMatchScoreFile"
#
# @title "A binary file holding match scores for each cell (probe/feature)"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \details{
#   Note that this class does \emph{not} assume a rectangular chip layout.
#   In other words, there is no concept of mapping a \emph{spatial}
#   location on the array to a cell index and vice versa.
#   The reason for this to be able to use this class also for
#   non-rectangular chip types.
# }
#
# @author "MR"

setConstructorS3("AromaCellMatchScoreFile", function(...) {
  extend(AromaCellTabularBinaryFile(...), "AromaCellMatchScoreFile");
})

setMethodS3("getFilenameExtension", "AromaCellMatchScoreFile", function(static, ...) {
  "acm";
}, static=TRUE)


setMethodS3("getExtensionPattern", "AromaCellMatchScoreFile", function(static, ...) {
  "[.](acm)$";
}, static=TRUE, protected=TRUE)


setMethodS3("getDefaultColumnNames", "AromaCellMatchScoreFile", function(this, ...) {
  c(sprintf("b%02d", seq(from=1, to=nbrOfColumns(this)-1)), "targetStrand");
})



setMethodS3("getDefaultExtension", "AromaCellMatchScoreFile", function(static, ...) {
  "cdf";
}, static=TRUE, protected=TRUE);


setMethodS3("byChipType", "AromaCellMatchScoreFile", function(static, chipType, tags=NULL, nbrOfCells=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'nbrOfCells':
  if (!is.null(nbrOfCells)) {
    nbrOfCells <- Arguments$getInteger(nbrOfCells, range=c(0,Inf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating ", class(static)[1])
  pathname <- findByChipType(static, chipType=chipType, tags=tags,
                                                 firstOnly=TRUE, ...);
  if (is.null(pathname)) {
      throw("Could not locate a file for this chip type: ",
                             paste(c(chipType, tags), collapse = ","));
  }
  verbose && cat(verbose, "Located file: ", pathname);
  res <- newInstance(static, pathname);
  verbose && print(verbose, res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validation?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(nbrOfCells)) {
    if (nbrOfCells(res) != nbrOfCells) {
      throw("The number of cells in the loaded ", class(static)[1], " does not match the expected number: ", nbrOfCells(res), " != ", nbrOfCells);
    }
  }

  verbose && exit(verbose);
  res;
}, static=TRUE)



setMethodS3("readMatchScores", "AromaCellMatchScoreFile", function(this, cells=NULL, drop=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    nbrOfCells <- length(cells);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Reading match scores");

  # Read data
  verbose && enter(verbose, "Reading data frame");
  res <- readDataFrame(this, rows=cells, columns=1, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # Flatten (data frame)
  dim <- dim(res);
  res <- unlist(res, use.names=FALSE);

  verbose && exit(verbose);

  res;
})


setMethodS3("updateMatchScores", "AromaCellMatchScoreFile", function(this, cells=NULL, scores, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, max=nbrOfCells);
    nbrOfCells <- length(cells);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Optimize
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove duplicated 'cells'
  keep <- which(!duplicated(cells));
  cells <- cells[keep];
  scores <- scores[keep];
  # Not needed anymore
  keep <- NULL;

  # Order by 'cells'
  srt <- sort(cells, method="quick", index.return=TRUE);
  o <- srt$ix;
  cells <- srt$x;
  # Not needed anymore
  srt <- NULL;
  scores <- scores[o];
  # Not needed anymore
  o <- NULL;

  lastColumn <- nbrOfColumns(this);
  this[cells,lastColumn] <- as.integer(scores);
})


setMethodS3("isMissing", "AromaCellMatchScoreFile", function(this, ...) {
  res <- readMatchScores(this, ..., positions=1, drop=TRUE);
  res <- (res == 0);
  res;
}, protected=TRUE)


setMethodS3("allocate", "AromaCellMatchScoreFile", function(static, ..., nbrOfCells, platform, chipType, footer=list()) {
  # Argument 'nbrOfCells':
  nbrOfCells <- Arguments$getInteger(nbrOfCells, range=c(1, 1000e6));

  # Argument 'platform':
  platform <- Arguments$getCharacter(platform, length=c(1,1));

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'footer':
  if (is.null(footer)) {
  } else if (!is.list(footer)) {
    throw("Argument 'footer' must be NULL or a list: ", class(footer)[1]);
  }

  footer <- c(
    list(
      createdOn=format(Sys.time(), "%Y%m%d %H:%M:%S", usetz=TRUE),
      platform=platform,
      chipType=chipType
    ),
    footer
  );

  nbrOfColumns <- 1L;
  NextMethod("allocate", nbrOfRows=nbrOfCells, types=rep("integer", times=1L), sizes=rep(1L, times=nbrOfColumns), footer=footer);
}, static=TRUE)


############################################################################
# HISTORY:
# 2010-02-18 [MR]
# o corrected the column names for AromaCellMatchscoreFile class
# 2009-02-10 [HB]
# o Added optional validation of number of cells to byChipType().
# o Static method byChipType() was not declared static.
# 2008-10-28 [MR]
# o Created from AromaCellSequenceFile.R.
############################################################################
