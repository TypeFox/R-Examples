###########################################################################/**
# @RdocClass "AromaCellCpgFile"
#
# @title "A binary file holding local CpG density for each cell (probe/feature)"
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
# @author "MR"
#*/###########################################################################
setConstructorS3("AromaCellCpgFile", function(...) {
  extend(AromaCellTabularBinaryFile(...), "AromaCellCpgFile");
})

setMethodS3("getFilenameExtension", "AromaCellCpgFile", function(static, ...) {
  "acc";
}, static=TRUE)


setMethodS3("getDefaultExtension", "AromaCellCpgFile", function(static, ...) {
  "acc";
}, static=TRUE)


setMethodS3("getDefaultColumnNames", "AromaCellCpgFile", function(this, ...) {
  c("cpgDensity");
}, protected=TRUE)


setMethodS3("readCpgs", "AromaCellCpgFile", function(this, cells=NULL, drop=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
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


setMethodS3("updateCpgs", "AromaCellCpgFile", function(this, cells=NULL, scores, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  nbrOfCells <- nbrOfCells(this);
  if (!is.null(cells)) {
    cells <- Arguments$getIndices(cells, range=c(1, nbrOfCells));
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


setMethodS3("isMissing", "AromaCellCpgFile", function(this, ...) {
  res <- readCpgs(this, ..., positions=1, drop=TRUE);
  res <- (res == 0);
  res;
}, protected=TRUE)


setMethodS3("allocate", "AromaCellCpgFile", function(static, ..., nbrOfCells, platform, chipType, footer=list()) {
  # Argument 'nbrOfCells':
  nbrOfCells <- Arguments$getInteger(nbrOfCells, range=c(1, 1000e6));

  # Argument 'platform':
  platform <- Arguments$getCharacter(platform);

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

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

  NextMethod("allocate", nbrOfRows=nbrOfCells, types="double", sizes=4L, footer=footer);
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2014-06-28
# o CLEANUP: Move byChipType() from classes AromaCellCpgFile and
#   AromaCellPositionFile to superclass AromaCellTabularBinaryFile.
# o Added argument 'nbrOfCells' to byChipType().
# o CLEANUP: Dropped non-used argument 'validate' from byChipType().
# 2010-02-19 [MR]
# o Added begin/end Rdoc comments so that they are compiled into Rd files.
# o Added file to SVN repository.
# o Modified header description, added semicolons.
# 2008-12-15 [MR]
# o Created from AromaCellPositionFile.R.
############################################################################
