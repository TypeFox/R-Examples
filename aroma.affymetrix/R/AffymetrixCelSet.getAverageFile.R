###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod getAverageFile
#
# @title "Calculates the mean and the standard deviation of the cell signal (intensity, standard deviation etc.) across the CEL set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{name}{The label of the calculated parameters.
#    If @NULL, a default name format \code{<prefix>-<mean>-<sd>} is used.}
#  \item{indices}{An @integer @vector specifying which cells to consider.
#    If \code{"remaining"}, only parameters for cells that have not been
#    are calculated.
#    If @NULL, all cells are used.}
#  \item{mean}{A @character of a @function specifying the function used
#    to calculate the average.}
#  \item{sd}{A @character of a @function specifying the function used
#    to calculate the standard deviation.}
#  \item{na.rm}{If @TRUE, @NAs are excluded before, otherwise not.}
#  \item{...}{Not used.}
#  \item{cellsPerChunk}{A @integer specifying the total number of cells
#    (across arrays) read into memory per chunk.}
#  \item{moreCells}{A @double scalar indicating if more or less cells
#    should be used per chunk.}
#  \item{force}{If @TRUE, parameters for cells already calculated are
#    recalculated, otherwise not.}
#  \item{verbose}{If @TRUE, progress details are printed, otherwise not.
#    May also be a @see "R.utils::Verbose" object.}
# }
#
# \value{
#   Returns an @see "AffymetrixCelSet" of the same class as the CEL set
#   averaged.
# }
#
# \details{
#   The parameter estimates are stored as a CEL file of the same class as
#   the data files in the set.  The CEL file is named \code{<name>.cel}
#   and placed in the directory of the set.
#   Currently there is no specific data class for this file, but the average
#   cell signals are stored as "intensities", the standard deviation of the
#   cell signals as "stddevs", and the number of data points used for each
#   estimate is stored as "pixels".
# }
#
# @author "HB, KS"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAverageFile", "AffymetrixCelSet", function(this, name=NULL, prefix="average", indices="remaining", field=c("intensities", "stdvs"), mean=c("median", "mean"), sd=c("mad", "sd"), na.rm=TRUE, g=NULL, h=NULL, ..., cellsPerChunk=moreCells*10^7/length(this), moreCells=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'field':
  field <- match.arg(field);

  # Argument 'mean':
  if (is.character(mean)) {
    mean <- match.arg(mean);
    meanName <- mean;
    if (mean == "mean") {
      mean <- base::rowMeans;
    } else if (mean == "median") {
      mean <- rowMedians;
    }
  } else if (is.function(mean)) {
    meanName <- "customMean";
  } else {
    throw("Argument 'mean' must be either a character or a function: ", mode(mean));
  }

  # Argument 'sd':
  if (is.character(sd)) {
    sd <- match.arg(sd);
    sdName <- sd;
    if (sd == "sd") {
      sd <- rowSds;
    } else if (sd == "mad") {
      sd <- rowMads;
    }
  } else if (is.function(sd)) {
    sdName <- "customSd";
  } else {
    throw("Argument 'sd' must be either a character or a function: ",
                                                           mode(sd));
  }

  # Argument 'name':
  if (is.null(name)) {
    key <- list(method="getAverageFile", class=class(this)[1],
                arrays=sort(getNames(this)), mean=meanName, sd=sdName);
    # assign mean and sd to an empty environment so that digest() doesn't
    # pick up any "promised" objects from the original environment.
    # A bit ad hoc, but it works for now. /2007-01-03
    key <- lapply(key, FUN=function(x) {
      if (is.function(x))
        environment(x) <- emptyenv();
      x;
    })
    id <- getChecksum(key);
    name <- sprintf("%s-%s-%s-%s,%s", prefix, field, meanName, sdName, id);
  }

  # Argument 'indices':
  df <- as.list(this)[[1]];
  nbrOfCells <- getHeader(df)$total;
  if (force) {
    if (identical(indices, "remaining")) {
      indices <- NULL;
    }
  }

  if (is.null(indices)) {
    indices <- 1:nbrOfCells;
  } else if (identical(indices, "remaining")) {
  } else {
    indices <- Arguments$getIndices(indices, max=nbrOfCells);
  }

  # Argument 'cellsPerChunk':
  cellsPerChunk <- Arguments$getInteger(cellsPerChunk, range=c(1,Inf));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving average cell signals across ", length(this), " arrays");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create CEL file to store the average array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a private filename (with a dot prefix) to make sure it is not
  # identified as a regular CEL file when the directory is scanned for files.
  filename <- sprintf(".%s.CEL", name);
  if (is.null(this$.averageFiles))
    this$.averageFiles <- list();
  res <- this$.averageFiles[[filename]];

  # Has file been deleted since last time?
  if (!is.null(res) && !isFile(res)) {
    warning("Will recalculate average file, because it seems to have been deleted since last time: ", getPathname(res));
    res <- NULL;
  }

  if (is.null(res)) {
    verbose && enter(verbose, "Obtaining an (existing or new) result file");

    # Searching for the output file in multiple directories
    path <- getPath(this);
    paths <- c(path);

    # Drop tags from root path?
    if (getOption(aromaSettings, "devel/dropRootPathTags", TRUE)) {
      path <- dropRootPathTags(path, depth=2, verbose=less(verbose, 5));
      paths <- c(paths, path);
      paths <- unique(paths);
    }

    verbose && cat(verbose, "Paths:");
    verbose && print(verbose, paths);
    verbose && cat(verbose, "Filename: ", filename);

    pathname <- NULL;
    for (kk in seq_along(paths)) {
      path <- paths[kk];
      verbose && enter(verbose, sprintf("Searching path #%d of %d", kk, length(paths)));

      verbose && cat(verbose, "Path: ", path);
      pathnameT <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
      verbose && cat(verbose, "Pathname: ", pathnameT);
      if (isFile(pathnameT)) {
        pathname <- pathnameT;
        verbose && cat(verbose, "Found an existing file.");
        verbose && exit(verbose);
        break;
      }

      verbose && exit(verbose);
    } # for (kk ...)
    verbose && cat(verbose, "Located pathname: ", pathname);


    if (isFile(pathname)) {
      verbose && enter(verbose, "Loading existing data file");
      verbose && cat(verbose, "Pathname: ", pathname);
      res <- newInstance(df, pathname);
      verbose && exit(verbose);
    } else {
      verbose && enter(verbose, "Creating CEL file to store average signals");
      path <- paths[length(paths)];

      verbose && cat(verbose, "Path: ", path);
      verbose && cat(verbose, "Filename: ", filename);
      pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
      verbose && cat(verbose, "Pathname: ", pathname);

      res <- createFrom(df, filename=pathname, path=NULL, methods="create", clear=TRUE, verbose=less(verbose));
      verbose && exit(verbose);
    } # if (isFile(pathname))

    verbose && exit(verbose);
    this$.averageFiles[[filename]] <- res;
  } # if (is.null(res))

  verbose && print(verbose, res);

  pathname <- getPathname(res);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify which indices to use
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(indices, "remaining")) {
    pixels <- .readCel(pathname, readIntensities=FALSE, readStdvs=FALSE,
                      readPixels=TRUE)$pixels;
    indices <- which(pixels == 0);
    # Not needed anymore
    pixels <- NULL; # Not needed anymore.
  }

  nbrOfIndices <- length(indices);

  # Nothing more to do?
  if (nbrOfIndices == 0)
    return(res);

  verbose && cat(verbose, "Number of cells to be updated: ", nbrOfIndices);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Estimate the mean and standard deviation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Since we might want to do this robustly, but also because we want to
  # estimate the standard deviation, for each cell we need all data across
  # arrays at once.  In order to this efficiently, we do this in chunks
  idxs <- 1:nbrOfIndices;
  head <- 1:cellsPerChunk;
  nbrOfChunks <- ceiling(nbrOfIndices / cellsPerChunk);
  verbose && cat(verbose, "Number cells per chunk: ", cellsPerChunk);

  # Get the pathnames of all CEL files to average
  pathnames <- lapply(this, FUN=getPathname);
  pathnames <- unlist(pathnames, use.names=FALSE);
  nbrOfArrays <- length(pathnames);

  if (!na.rm)
    n <- rep(nbrOfArrays, length=cellsPerChunk);
  count <- 1;
  while (length(idxs) > 0) {
    verbose && enter(verbose, "Fitting chunk #", count, " of ", nbrOfChunks);
    if (length(idxs) < cellsPerChunk) {
      head <- 1:length(idxs);
      if (!na.rm)
        n <- rep(nbrOfArrays, length=length(idxs));
    }

    # The indices to be used in this chunk
    ii <- idxs[head];
    verbose && cat(verbose, "Chunk size: ", length(ii));

    verbose && enter(verbose, "Reading data");
#    X <- readCelIntensities(pathnames, indices=indices[ii]);
    readIntensities <- (field == "intensities");
    readStdvs <- (field == "stdvs");
    # TODO: Ideally, affxparser::readCel() should support
    # multiple filenames turning every data fields into a
    # matrix. /HB 2007-01-07
    X <- matrix(as.double(NA), nrow=length(ii), ncol=nbrOfArrays);
    for (kk in seq_len(nbrOfArrays)) {
      X[,kk] <- .readCel(filename = pathnames[kk],
                        indices = indices[ii],
                        readIntensities = readIntensities,
                        readHeader = FALSE,
                        readStdvs = readStdvs,
                        readPixels = FALSE,
                        readXY = FALSE,
                        readOutliers = FALSE,
                        readMasked = FALSE,
                        ...,
                        verbose = (verbose - 1))[[field]];
    }
    verbose && exit(verbose);

    if (!is.null(g)) {
      verbose && enter(verbose, "Transforming data using y = g(x)");
      X <- g(X);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Estimating averages and standard deviations");
    if (na.rm)
      n <- base::apply(X, MARGIN=1, FUN=function(x) { sum(!is.na(x)) });

    # Calculate the mean signal
    mu <- mean(X, na.rm=na.rm);          # Special mean()!
    # Calculate the standard deviation of the signals
    sigma <- sd(X, mean=mu, na.rm=na.rm);   # Special sd()!

    verbose && exit(verbose);

    if (!is.null(h)) {
      verbose && enter(verbose, "Back-transforming estimates using x = h(y)");
      mu <- h(mu);
      sigma <- h(sigma);
      verbose && exit(verbose);
    }

    # Write estimates to result file
    verbose && enter(verbose, "Writing estimates");
    .updateCel(pathname, indices=indices[ii], intensities=mu, stdvs=sigma, pixels=n);
    verbose && exit(verbose);

    # Not needed anymore
    mu <- sigma <- NULL;

    # Next chunk...
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    gc();
    verbose && exit(verbose);
  } # while()

  verbose && exit(verbose);

  res;
})



############################################################################
# HISTORY:
# 2011-02-28
# o Now getAverageFile() first tries to locate an existing result file
#   in multiple root paths.  If not found, it creates a new one.
# 2011-02-24
# o GENERALIZATION: Now getAverageFile() for AffymetrixCelSet drops tags
#   from the output root path (if 'devel/dropRootPathTags' setting is TRUE).
# 2010-07-21
# o CLEAN UP: Now getAverageFile() for AffymetrixCelSet no longer writes
#   debug information to ${Rcache}/aroma.affymetrix/idChecks/.
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: getData(), getAverageFile().
# 2007-09-25
# o Now getAverageFile() will detect if an average file has been deleted
#   between calls and recalculate it.
# 2007-03-16
# o BUG FIX: getAverageFile() of AffymetrixCelSet would average the wrong
#   set of cells if argument 'indices' was different from NULL.
# 2007-01-07
# o BUG FIX: In KS's update of getAverageFile() to support averaging
#   over other fields than intensities, argument 'indices' was missing
#   in the readCel() call making the function fail when processed chunk
#   by chunk. /HB
# 2006-11-07
# o Now getAverageFile() uses rowMedians() of R.native if available,
#   otherwise a local version utilizing apply(). Same for rowMads().
# 2006-10-24
# o Added getAverageLog() and getAverageAsinh().
# o Added transforms and anti-transforms g() and h() to getAverageFile().
# o Changed the defaults from mean to median, and sd to mad for
#   getAverageFile().
# o Added Rdoc comments to getAverageFile().
# 2006-08-27
# o Added getAverageFile().
############################################################################
