###########################################################################/**
# @RdocClass FileMatrix
# @alias FileByteMatrix
# @alias FileShortMatrix
# @alias FileIntegerMatrix
# @alias FileFloatMatrix
# @alias FileDoubleMatrix
#
# @title "Class representing a persistent matrix stored in a file"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "AbstractFileArray".}
#  \item{nrow, ncol}{The number of rows and columns of the matrix.}
#  \item{rownames, colnames}{Optional row and column names.}
#  \item{byrow}{If @TRUE, data are stored row by row, otherwise column
#      by column.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
#
# }
#
# \details{
#   The purpose of this class is to be able to work with large matrices
#   in \R without being limited by the amount of memory available.
#   Matrices are kept on the file system and elements are read and written
#   whenever queried.  The purpose of the class is \emph{not} to provide
#   methods for full matrix operations, but instead to be able to work
#   with subsets of the matrix at each time.
#
#   For more details, @see "AbstractFileArray".
# }
#
# \section{Column by column or row by row?}{
#   If the matrix elements are to be accessed more often along rows,
#   store data row by row, otherwise column by column.
# }
#
# \section{Supported data types}{
#   The following subclasses implement support for various data types:
#   \itemize{
#    \item \code{FileByteMatrix} (1 byte per element),
#    \item \code{FileShortMatrix} (2 bytes per element),
#    \item \code{FileIntegerMatrix} (4 bytes per element),
#    \item \code{FileFloatMatrix} (4 bytes per element), and
#    \item \code{FileDoubleMatrix} (8 bytes per element).
#   }
# }
#
# @examples "../incl/FileMatrix.Rex"
#
# @author
#
# @visibility public
#*/###########################################################################
setConstructorS3("FileMatrix", function(..., nrow=NULL, ncol=NULL, rownames=NULL, colnames=NULL, byrow=FALSE) {

  dim <- c(nrow, ncol);
  dimnames <- list(rownames, colnames);
  dimOrder <- 1:2;
  if (byrow)
    dimOrder <- 2:1;
  extend(AbstractFileArray(..., dim=dim, dimnames=dimnames, dimOrder=dimOrder), "FileMatrix");
})



###########################################################################/**
# @RdocMethod as.character
#
# @title "Returns a short string describing the file matrix"
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
#  Returns a @character string.
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
setMethodS3("as.character", "FileMatrix", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf(" byrow=%s.", getByRow(this)));
  s <- paste(s, collapse="");

  s;
})


###########################################################################/**
# @RdocMethod getByRow
#
# @title "Checks if elements are stored row by row or not"
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
#  Returns @TRUE if the elements are stored row by row, otherwise @FALSE.
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
setMethodS3("getByRow", "FileMatrix", function(this, ...) {
  as.logical(diff(this$header$dimOrder) < 0);
})



###########################################################################/**
# @RdocMethod nrow
#
# @title "Gets the number of rows of the matrix"
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
#  Returns a @double.
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
setMethodS3("nrow", "FileMatrix", function(this, ...) {
  dim(this)[1];
})



###########################################################################/**
# @RdocMethod ncol
#
# @title "Gets the number of columns of the matrix"
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
#  Returns a @double.
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
setMethodS3("ncol", "FileMatrix", function(this, ...) {
  dim(this)[2];
})


###########################################################################/**
# @RdocMethod rownames
#
# @title "Gets the row names of a file matrix"
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
#  Returns a @character @vector, or @NULL.
# }
#
# @author
#
# \seealso{
#   @seemethod "colnames".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("rownames", "FileMatrix", function(this, ...) {
  dimnames(this)[[1]];
})


###########################################################################/**
# @RdocMethod colnames
#
# @title "Gets the column names of a file matrix"
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
#  Returns a @character @vector, or @NULL.
# }
#
# @author
#
# \seealso{
#   @seemethod "rownames".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("colnames", "FileMatrix", function(this, ...) {
  dimnames(this)[[2]];
})



setMethodS3("getColumnOffset", "FileMatrix", function(this, cols, offset=getDataOffset(this), bytesPerCell=getBytesPerCell(this), ...) {
  if (getByRow(this)) {
    offset <- offset + bytesPerCell*(cols-1);
  } else {
    offset <- offset + bytesPerCell*nrow(this)*(cols-1);
  }
  offset;
}, protected=TRUE)

setMethodS3("getRowOffset", "FileMatrix", function(this, rows, offset=getDataOffset(this), bytesPerCell=getBytesPerCell(this), ...) {
  if (getByRow(this)) {
    offset <- offset + bytesPerCell*ncol(this)*(rows-1);
  } else {
    offset <- offset + bytesPerCell*(rows-1);
  }
  offset;
}, protected=TRUE)



setMethodS3("getOffset", "FileMatrix", function(this, rows, cols, offset=getDataOffset(this), bytesPerCell=getBytesPerCell(this), ...) {
  if (length(rows) > length(cols)) {
    cols <- rep(cols, length.out=length(rows));
  } else if (length(rows) < length(cols)) {
    rows <- rep(rows, length.out=length(cols));
  }

  if (getByRow(this)) {
    idx <- ncol(this)*(rows-1) + (cols-1);
  } else {
    idx <- nrow(this)*(cols-1) + (rows-1);
  }

  offset <- offset + bytesPerCell*idx;

  offset;
}, protected=TRUE)



setMethodS3("writeValues", "FileMatrix", function(this, rows, cols, values, order=FALSE, ...) {
  con <- this$con;

  # Transform values to correct storage mode
  storage.mode(values) <- getStorageMode(this);

  # Get offsets for each of the values
  offsets <- getOffset(this, rows=rows, cols=cols);

  if (length(offsets) != length(values)) {
    throw("Number of indices and values does not match.");
  }

  # Order offsets to optimize caching?
  if (order) {
    o <- order(offsets);
    offsets <- offsets[o];
    values <- values[o];
  }

  what <- values[1];
  size <- getBytesPerCell(this);

  ints <- seqToIntervals(offsets);
  valueOffset <- 0;

  for (rr in seq(length=nrow(ints))) {
    from <- ints[rr,1];
    to <- ints[rr,2];
    n <- to-from+1;
    idx <- valueOffset + 1:n;
    tmp <- values[idx];
    valueOffset <- valueOffset + n;
    seek(con=con, where=from, rw="write");
    writeBin(con=con, tmp, size=size);
  }
}, protected=TRUE)


setMethodS3("readValues", "FileMatrix", function(this, rows, cols, order=FALSE, ...) {
  con <- this$con;

  # Get offsets for each of the values
  offsets <- getOffset(this, rows=rows, cols=cols);

  # Allocate vector for values;
  values <- vector(mode=getStorageMode(this), length=length(offsets));

  if (order) {
    o <- order(offsets);
    offsets <- offsets[o];
  }

  what <- vector(mode=getStorageMode(this), length=1);
  size <- getBytesPerCell(this);

  ints <- seqToIntervals(offsets);
  valueOffset <- 0;
  for (rr in seq(length=nrow(ints))) {
    from <- ints[rr,1];
    to <- ints[rr,2];
    n <- to-from+1;
    idx <- valueOffset + 1:n;
    valueOffset <- valueOffset + n;
    seek(con=con, where=from, rw="read");
    tmp <- readBin(con=con, what=what, n=n, size=size);
    values[idx] <- tmp;
  }

  values;
}, protected=TRUE)


setMethodS3("readFullMatrix", "FileMatrix", function(this, ...) {
  values <- readAllValues(this);

  if (getByRow(this)) {
    dim(values) <- rev(dim(this));
    names <- rownames(this);
    if (length(names) > 0)
      colnames(values) <- names;
    names <- colnames(this);
    if (length(names) > 0)
      rownames(values) <- names;
  } else {
    dim(values) <- dim(this);
    names <- rownames(this);
    if (length(names) > 0)
      rownames(values) <- names;
    names <- colnames(this);
    if (length(names) > 0)
      colnames(values) <- names;
  }

  values;
}, protected=TRUE)


setMethodS3("getMatrixIndicies", "FileMatrix", function(this, i, j, ...) {
  nrow <- this$header$dim[1];
  ncol <- this$header$dim[2];

  if (missing(i)) {
    i <- 1:nrow;
  } else if (any(i < 1 | i > nrow)) {
    throw("Argument 'i' is out of range.");
  }

  if (missing(j)) {
    j <- 1:ncol;
  } else if (any(j < 1 | j > ncol)) {
    throw("Argument 'j' is out of range.");
  }

  ni <- length(i);
  nj <- length(j);

  dataOffset <- getDataOffset(this);
  bytesPerCell <- getBytesPerCell(this);

  if (getByRow(this)) {
    # Matrix elements stored row by row
    iOffsets <- bytesPerCell*ncol*(i - 1);
    jOffsets <- bytesPerCell*(j - 1);
    iOffsets <- dataOffset + iOffsets;
  } else {
    # Matrix elements stored column by column
    iOffsets <- bytesPerCell*(i - 1);
    jOffsets <- bytesPerCell*nrow*(j - 1);
    jOffsets <- dataOffset + jOffsets;
  }

  list(
    nrow = nrow,
    ncol = ncol,
    byrow = getByRow(this),
    i = i,
    j = j,
    ni = ni,
    nj = nj,
    bytesPerCell = bytesPerCell,
    iOffsets = iOffsets,
    jOffsets = jOffsets
  )
}, protected=TRUE)



###########################################################################/**
# @RdocMethod "["
#
# @title "Gets a subset of elements of a file matrix"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{A @numeric @vector or row indices.}
#   \item{j}{A @numeric @vector or column indices.}
#   \item{drop}{If @TRUE, if only a single row or column is retrived,
#     it is returned as a @vector.}
# }
#
# \value{
#  Returns a @matrix (or possibly a @vector).
# }
#
# @author
#
# \seealso{
#   @seemethod "[<-".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("[", "FileMatrix", function(this, i, j, drop=FALSE) {
  # Calculate row and column index offsets
  res <- getMatrixIndicies(this, i, j);
  i <- res$i;
  j <- res$j;
  ni <- res$ni;
  nj <- res$nj;
  ncol <- res$ncol;
  nrow <- res$nrow;
  byrow <- res$byrow;
  size <- res$bytesPerCell;
  res <- NULL; # Not needed anymore

  if (byrow) {
    i <- (i-as.integer(1))*ncol;
  } else {
    j <- (j-as.integer(1))*nrow;
  }

  # Get the location of all cells to be read
  idxs <- outer(i, j, FUN="+");
  idxs <- as.vector(idxs);
  i <- j <- NULL; # Not needed anymore

  # Move to the data section
  con <- this$con;
  seek(con=con, where=getDataOffset(this), rw="read");

  mode <- getStorageMode(this);
  res <- readBinFragments(con=con, what=mode, size=size,
                          idxs=idxs, origin="current");
  dim(res) <- c(ni, nj);

  # Dimension names?
  if (!drop) {
    names <- rownames(this);
    if (length(names) > 0)
      rownames(res) <- names;
    names <- colnames(this);
    if (length(names) > 0)
      colnames(res) <- names;
  }

  # Drop indices?
  if (drop) {
    if (ni <= 1 || nj <= 1)
      dim(res) <- NULL;
  }

  res;
}) # "["()


###########################################################################/**
# @RdocMethod as.matrix
#
# @title "Returns the elements of a file matrix as an R matrix"
#
# \description{
#  @get "title", that is, imported into memory (if possible).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @matrix.
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
setMethodS3("as.matrix", "FileMatrix", function(x, ...) {
  # To please R CMD check
  this <- x;

  this[];
})


###########################################################################/**
# @RdocMethod "[<-"
#
# @title "Assigns values to a subset of elements of a file matrix"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{i}{A @numeric @vector or row indices.}
#   \item{j}{A @numeric @vector or column indices.}
#   \item{value}{Values to be assigned to the selected elements.}
# }
#
# \value{
#  Returns itself.
# }
#
# @author
#
# \seealso{
#   @seemethod "[".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("[<-", "FileMatrix", function(this, i, j, value) {
  # Calculate row and column index offsets
  res <- getMatrixIndicies(this, i, j);
  i <- res$i;
  j <- res$j;
  ncol <- res$ncol;
  nrow <- res$nrow;
  byrow <- res$byrow;
  size <- res$bytesPerCell;
  res <- NULL; # Not needed anymore

  if (byrow) {
    i <- (i-as.integer(1))*ncol;
  } else {
    j <- (j-as.integer(1))*nrow;
  }

  # Get the location of all cells to be read
  idxs <- outer(i, j, FUN="+");
  idxs <- as.vector(idxs);
  i <- j <- NULL; # Not needed anymore

  # Coerce input data
  storage.mode(value) <- getStorageMode(this);

  # Move to the data section
  con <- this$con;
  seek(con=con, where=getDataOffset(this), rw="write");

  # Write the data
  nValue <- length(value);
  nIdxs <- length(idxs);
  if (nValue < nIdxs) {
    value <- rep(value, length.out=nIdxs);
  } else if (nValue > nIdxs) {
    value <- value[seq(length=nIdxs)];
  }

  writeBinFragments(con=con, value, size=size, idxs=idxs);

  this;
})



###########################################################################/**
# @RdocMethod rowSums
#
# @title "Calculates the sum for each row"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{na.rm}{If @TRUE, @NA values are excluded.}
#   \item{doCount}{If @TRUE, the number of included values are counted
#      and returned as an attribute.}
#   \item{rows}{An @integer @vector of rows for which the sum should
#      be calculated.  If @NULL, all rows are considered.}
#   \item{columns}{An @integer @vector of columns for which the sum should
#      be calculated.  If @NULL, all columns are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "rowMeans".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("rowSums", "FileMatrix", function(x, na.rm=FALSE, doCount=FALSE, rows=NULL, columns=NULL, ...) {
  # To please R CMD check
  this <- x;

  nrow <- nrow(this);
  ncol <- ncol(this);

  if (is.null(rows)) {
    rows <- seq(length=nrow);
  } else {
    nrow <- length(rows);
  }

  if (is.null(columns)) {
    columns <- seq(length=ncol);
  } else {
    ncol <- length(columns);
  }

  if (getByRow(this)) {
    sums <- rep(NA, length=nrow);
    counts <- rep(0, length=nrow);
    if (na.rm) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Remove NAs before summing
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      for (rr in rows) {
        values <- this[rr,columns,drop=TRUE];
        ok <- !is.na(values);
        if (doCount)
          counts[rr] <- sum(ok);
        sums[rr] <- sum(values[ok], na.rm=FALSE);
        values <- NULL; # Not needed anymore
      } # for (rr in ...)
    } else {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Sum with NAs
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      for (rr in rows) {
        sums[rr] <- sum(this[rr,columns,drop=TRUE], na.rm=FALSE);
      } # for (rr in ...)
    } # if (na.rm)
  } else {
    sums <- NULL;
    if (na.rm) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Remove NAs before summing
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      for (cc in columns) {
        values <- this[rows,cc,drop=TRUE];
        if (cc == columns[1]) {
          sums <- values;
          currNAs <- is.na(sums);
          if (doCount)
            counts <- as.integer(!currNAs);
        } else {
          # NAs in current
          newNAs <- is.na(values);

          # (i) Sum elements with NAs can be assign with the new values
          sums[currNAs] <- values[currNAs];

          # NAs after adding this column
          currNAs <- (currNAs & newNAs);

          # (ii) Remaining non-NAs elements can be added
          sums[!currNAs] <- sums[!currNAs] + values[!currNAs];

          if (doCount)
            counts <- counts + as.integer(!newNAs);
        }
        values <- NULL; # Not needed anymore
      } # for (cc in ...)
    } else {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Sum with NAs
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      for (cc in columns) {
        values <- this[rows,cc,drop=TRUE];
        if (cc == columns[1]) {
          sums <- values;
        } else {
          sums <- sums + values;
        }
        values <- NULL; # Not needed anymore
      } # for (cc in ...)
    } # if (na.rm)
  } # if (getByRow(this))

  # Return counts too?
  if (doCount) {
    if (!na.rm)
      counts <- rep(ncol, length=nrow);
    attr(sums, "counts") <- counts;
  }

  sums;
})



###########################################################################/**
# @RdocMethod rowMeans
#
# @title "Calculates the means for each row"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @seemethod "rowSums".}
#   \item{doCount}{If @TRUE, the number of included values are counted
#      and returned as an attribute.}
# }
#
# \value{
#  Returns a @numeric @vector.
# }
#
# @author
#
# \seealso{
#   @seemethod "rowSums".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("rowMeans", "FileMatrix", function(x, ..., doCount=FALSE) {
  # To please R CMD check
  this <- x;

  res <- rowSums(this, doCount=TRUE, ...);
  res <- res/attr(res, "counts");
  if (!doCount)
    attr(res, "counts") <- NULL;
  res;
})


############################################################################
# HISTORY:
# 2010-11-08
# o Now "["() for FileMatrix explicitly specifies origin="current" in the
#   call to readBinFragment(), which is an argument added in R.utils 1.5.7.
# 2007-08-22
# o Added special as.character() for FileMatrix.
# o BUG FIX: "["() and "[<-"() were totally broken.  Now they utilize the
#   new readBinFragments() and writeBinFragments(), respectively.
#   Eventually this should be implement in AbstractFileArray as a generic
#   solution.
# 2006-05-09
# o Added Rdoc comments.
# 2006-02-27
# o Now inheriting from class AbstractFileArray.  Many of the low-level
#   functions have been moved to that class.
# 2006-02-18
# o Now all indices are reprensented as doubles instead of integers.
#   This make the dimension and/or number of elements in the matrix limited
#   only by the file system.
# o Added some Rdoc comments.
# 2006-01-23
# o First benchmarking:
#    a) Reading 90 Affymetrix Xba SNP chips, each with 2,500,000 probes,
#       takes approximately 15.0 minutes.
#    b) Then quantile normalizing (returning a copy) takes approximately
#       15.0 minutes too.
#    c) Plotting density functions for these takes 4.0 minutes.
#   (a)+(b) takes exactly 30 minutes.
# o Added rowSums() and rowMeans().
# o Now clone() identifies the clone number, and searches for the first
#   non-used clone number above that.
# o Added some simple Rdoc comments with example code.
# o TO DO: If 'j' is not ordered for X[i,j] <- ... with byrow=TRUE, wrong
#   columns are assigned.
# o Most works now.
# 2006-01-22
# o Created.
############################################################################
