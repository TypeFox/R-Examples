###########################################################################/**
# @RdocClass "AromaCellPositionFile"
#
# @title "A binary file holding chromosome/position for each cell"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to constructor of
#             @see "AromaCellTabularBinaryFile".}
# }
#
# \details{
#   Note that this class does \emph{not} assume a rectangular chip layout.
#   In other words, there is no concept of mapping a \emph{spatial}
#   location on the array to a cell index and vice versa.
#   The reason for this to be able to use this class also for
#   non-rectangular chip types.
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaCellPositionFile", function(...) {
  extend(AromaCellTabularBinaryFile(...), "AromaCellPositionFile");
})

setMethodS3("getFilenameExtension", "AromaCellPositionFile", function(static, ...) {
  "acp";
}, static=TRUE)

setMethodS3("getDefaultExtension", "AromaCellPositionFile", function(static, ...) {
  "acp";
}, static=TRUE)


setMethodS3("getExtensionPattern", "AromaCellPositionFile", function(static, ...) {
  "[.](acp)$";
}, static=TRUE, protected=TRUE)


setMethodS3("getDefaultColumnNames", "AromaCellPositionFile", function(this, ...) {
  c("chromosome", "position");
}, protected=TRUE)


setMethodS3("readPositions", "AromaCellPositionFile", function(this, cells=NULL, drop=FALSE, ..., verbose=FALSE) {
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
  res <- readDataFrame(this, rows=cells, columns=2, verbose=less(verbose, 5));
  verbose && exit(verbose);

  # Flatten (data frame)
  dim <- dim(res);
  res <- unlist(res, use.names=FALSE);

  verbose && exit(verbose);

  res;
})


setMethodS3("updatePositions", "AromaCellPositionFile", function(this, cells=NULL, scores, ..., verbose=FALSE) {
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


setMethodS3("isMissing", "AromaCellPositionFile", function(this, ...) {
  res <- readPositions(this, ..., positions=1, drop=TRUE);
  res <- (res == 0);
  res;
}, protected=TRUE)


setMethodS3("allocate", "AromaCellPositionFile", function(static, ..., nbrOfCells, platform, chipType, footer=list()) {
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

  NextMethod("allocate", nbrOfRows=nbrOfCells, types=rep("integer", times=2L), sizes=c(1L,4L), footer=footer);
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2014-06-28
# o CLEANUP: Move byChipType() from classes AromaCellCpgFile and
#   AromaCellPositionFile to superclass AromaCellTabularBinaryFile.
# 2009-02-16 [HB]
# Removed argument 'validate' from byChipType() of AromaCellPositionFile.
# 2009-02-10 [HB]
# o Added optional validation of number of cells to byChipType().
# o Static method byChipType() was not declared static.
# 2008-12-09 [MR]
# o Created from AromaCellMatchScoresFile.R.
############################################################################
