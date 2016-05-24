###########################################################################/**
# @RdocClass GenericTabularFile
#
# @title "The abstract GenericTabularFile class"
#
# \description{
#  @classhierarchy
#
#  A TabularTextFile is an object refering to a tabular text file
#  on a file system containing data in a tabular format.
#  Methods for reading all or a subset of the tabular data exist.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenericDataFile".}
#   \item{.verify, verbose}{(Internal only) If @TRUE, the file is
#      verified while the object is instantiated by the constructor.
#      The verbose argument is passed to the verifier function.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
#   An object of this class is typically part of an
#   @see "GenericTabularFileSet".
# }
#*/###########################################################################
setConstructorS3("GenericTabularFile", function(..., .verify=TRUE, verbose=FALSE) {
  this <- extend(GenericDataFile(...), c("GenericTabularFile", uses("ColumnNamesInterface")));

  if (.verify) {
    verify(this, ..., verbose=verbose);
  }

  this;
}, abstract=TRUE)


setMethodS3("as.character", "GenericTabularFile", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  s <- c(s, sprintf("Number of data rows: %d", nbrOfRows(this, fast=TRUE)));
  s;
}, protected=TRUE)



setMethodS3("verify", "GenericTabularFile", function(this, ..., verbose=FALSE) {
  # Nothing to do?
  pathname <- getPathname(this);
  if (is.null(pathname) || is.na(pathname))
    return(invisible(this));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Validating file contents");

  tryCatch({
    data <- readDataFrame(this, rows=1:10, verbose=verbose);
  }, error = function(ex) {
    throw("File format error of the tabular file ('", getPathname(this), "'): ", ex$message);
  })

  verbose && exit(verbose);

  invisible(this);
}, private=TRUE)





###########################################################################/**
# @RdocMethod nbrOfRows
# @alias nbrOfColumns.GenericTabularFile
#
# @title "Gets the number of data rows"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer.
# }
# @author
#
# \seealso{
#   @seemethod "dim".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("nbrOfRows", "GenericTabularFile", abstract=TRUE);



setMethodS3("nbrOfColumns", "GenericTabularFile", function(this, ...) {
  ncols <- NextMethod("nbrOfColumns");
  if (!is.na(ncols)) return(ncols);
  data <- readDataFrame(this, colClasses=NULL, rows=1L);
  ncol(data);
})



###########################################################################/**
# @RdocMethod dim
#
# @title "Gets the dimension of data table"
#
# \description{
#  @get "title", which is the number of rows and the number of columns.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns an @integer @vector of length two.
# }
# @author
#
# \seealso{
#   @seemethod "nbrOfRows".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("dim", "GenericTabularFile", function(x) {
  # To please R CMD check.
  this <- x;

  c(nbrOfRows(this), nbrOfColumns(this));
}, appendVarArgs=FALSE)






###########################################################################/**
# @RdocMethod readDataFrame
#
# @title "Reads the tabular data as a data frame"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readDataFrame", "GenericTabularFile", abstract=TRUE);



###########################################################################/**
# @RdocMethod readColumns
#
# @title "Reads a subset of the columns as a data frame"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns a @data.frame.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("readColumns", "GenericTabularFile", abstract=TRUE);




###########################################################################/**
# @RdocMethod extractMatrix
#
# @title "Reads one of the columns"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{column}{An @integer specifying the column to read.}
#   \item{drop}{If @TRUE, a @vector is returned,
#     otherwise a one-column @matrix.}
#   \item{...}{Additional arguments passed to @seemethod "readColumns".}
#   \item{verbose}{A @logical or a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns a Jx1 @matrix, or if \code{drop=TRUE} a @vector of length J.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("extractMatrix", "GenericTabularFile", function(this, column=1L, drop=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfColumns <- nbrOfColumns(this);

  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting data as a single-column matrix");

  # Read data as data frame
  data <- readColumns(this, columns=column, ..., verbose=less(verbose, 5));
  # Drop dimension
  data <- data[,1,drop=TRUE];

  verbose && cat(verbose, "Raw data frame read:");
  verbose && str(verbose, data);

  # Coerce into a matrix?
  if (!drop) {
    data <- as.matrix(data);
    colnames(data) <- getName(this);
  } else {
    verbose && cat(verbose, "Dropping singleton dimensions");
  }

  verbose && cat(verbose, "Result:");
  verbose && str(verbose, data);

  verbose && exit(verbose);

  data;
})


setMethodS3("[", "GenericTabularFile", function(this, i=NULL, j=NULL, drop=FALSE) {
  # Argument 'drop':
  drop <- Arguments$getLogical(drop);

  # Read data
  if (missing(j) || is.null(j)) {
    data <- readColumns(this, rows=i);
  } else {
    data <- readColumns(this, rows=i, columns=j);
  }

  # Drop dimensions?
  if (drop) {
    if (ncol(data) == 1L) {
      data <- data[,1L];
    } else if (nrow(data) == 1L) {
      data <- data[1L,];
    }
  }

  data;
}, protected=TRUE)


setMethodS3("head", "GenericTabularFile", function(x, n=6L, ...) {
  stopifnot(length(n) == 1L);
  nrow <- nrow(x);
  if (n < 0L) {
    n <- max(nrow + n, 0L);
  } else {
    n <- min(n, nrow);
  }
  rows <- seq_len(n);
  x[rows,, drop=FALSE];
})


setMethodS3("tail", "GenericTabularFile", function(x, n=6L, ...) {
  stopifnot(length(n) == 1L);
  nrow <- nrow(x);
  if (n < 0L) {
    n <- max(nrow + n, 0L);
  } else {
    n <- min(n, nrow);
  }
  rows <- seq.int(to=nrow, length.out=n);
  x[rows,, drop=FALSE];
})



############################################################################
# HISTORY:
# 2013-12-18
# o Added nbrOfColumns() for GenericTabularFile, which, if the number
#   of columns cannot be inferred from the column names, will fall back
#   to read the first row of data and use that as the number of columns.
# 2012-12-08
# o GENERALIZATION: Moved "["() to GenericTabularFile (from
#   TabularTextFile) and made it utilize readColumns().
# 2012-11-02
# o CLEANUP: Dropped all methods that are now in ColumnNamesInterface, e.g.
#   getColumnNames(), setColumnNames(), getColumnNamesTranslator().
# o Now GenericTabularFile implements ColumnNamesInterface.
# o CLEANUP: Moved nbrOfColumns() to ColumnNamesInterface.
# 2012-11-01
# o Added setColumnNames() for GenericTabularFile, which utilizes
#   setColumnNamesTranslator().
# 2012-09-01
# o CONSISTENCY: Now extractMatrix() for GenericTabularFile adds column
#   names just as ditto for GenericTabularFileSet does.
# 2010-08-16
# o Added some Rdoc comments.
# 2009-10-30
# o ROBUSTIFICATION: Now translateColumnNames() of GenericTabularFile throws
#   an exception if some fullnames were translated into NA.
# 2008-05-12
# o Added extractMatrix().
# o BUG FIX: getReadArguments() did not infer column classes if there was
#   no header to read but the column names was manually set.
# o BUG FIX: readDataFrame() did not read the first data row if there was
#   no column header; it was eaten up by a preceeding readHeader().
# 2008-04-29
# o Added readLines(), nbrOfLines(), nbrOfRows() and dim().
# o Now readDataFrame() keeps the row names if arguments rows != NULL.
# 2008-04-25
# o Now argument 'verbose' of the constructor is passed to verfity().
# 2008-04-24
# o Added argument 'rows' to readDataFrame() for GenericTabularFile.
# 2008-04-14
# o Renamed readData() to readDataFrame() for GenericTabularFile.
# 2008-03-22
# o Added {get|set}ColumnNameTranslator().
# 2008-03-18
# o Now any '...' arguments to getReadArguments() override the inferred
#   read arguments, e.g. na.strings="NA".
# 2008-02-27
# o Since 'affy' defines standardGeneric("colnames") and because S3 methods
#   are not found by such S4 generic functions, we avoid that method name,
#   and instead use getColumnNames().
# 2007-09-16
# o Removed all 'camelCaseNames' arguments.  Now column names are decided
#   by getColumnNames() and translateColumnNames(), which can be overridden.
# 2007-09-14
# o Extracted from AffymetrixTabularFile.
# 2007-09-10
# o Created from AffymetrixCsvGenomeInformation.R.
############################################################################
