#########################################################################/**
# @RdocDefault readApdUnits
#
# @title "Reads Affymetrix probe data (APD) as units (probesets)"
#
# @synopsis
#
# \description{
#   @get "title" by using the unit and group definitions in the
#   corresponding Affymetrix CDF file.
#
#   If more than one APD file is read, all files are assumed to be of
#   the same chip type, and have the same read map, if any.
#   It is not possible to read APD files of different types at the
#   same time.
# }
#
# \arguments{
#   \item{filenames}{The filenames of the APD files.  All APD files must
#     be of the same chip type.}
#   \item{units}{An @integer @vector of unit indices specifying which units
#     to be read.  If @NULL, all units are read.}
#   \item{...}{Additional arguments passed to @see "readApd".}
#   \item{transforms}{A @list of exactly \code{length(filenames)}
#     @functions.  If @NULL, no transformation is performed.
#     Values read are passed through the corresponding transform
#     function before being returned.}
#   \item{cdf}{A @character filename of a CDF file, or a CDF @list
#     structure.  If @NULL, the CDF file is searched for by
#     @see "affxparser::findCdf" first starting from the current directory
#     and then from the directory where the first APD file is.}
#   \item{stratifyBy}{Argument passed to low-level method
#     @see "affxparser::readCdfCellIndices".}
#   \item{addDimnames}{If @TRUE, dimension names are added to arrays,
#     otherwise not.  The size of the returned APD structure in bytes
#     increases by 30-40\% with dimension names.}
#   \item{readMap}{A @vector remapping cell indices to file indices.
#     If \code{"byMapType"}, the read map of type according to APD header
#     will be search for and read.  It is much faster to specify the
#     read map explicitly compared with searching for it each time.
#     If @NULL, no map is used.}
#   \item{dropArrayDim}{If @TRUE and only one array is read, the elements of
#     the group field do \emph{not} have an array dimension.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   A named @list where the names corresponds to the names of the units
#   read.  Each element of the @list is in turn a @list structure
#   with groups (aka blocks).
# }
#
# \section{Speed}{
#   Since the cell indices are semi-randomized across the array and
#   with units (probesets), it is very unlikely that the read will
#   consist of subsequent cells (which would be faster to read).
#   However, the speed of this method, which uses @see "R.huge::FileVector"
#   to read data, is comparable to the speed of
#   @see "affxparser::readCelUnits", which uses the Fusion SDK
#   (@see "affxparser::readCel") to read data.
# }
#
# @author
#
# @examples "../incl/readApdUnits.Rex"
#
# \seealso{
#   To read CEL units, @see "affxparser::readCelUnits".
#   Internally, the @see "readApd" method is used for read probe data,
#   and @see "readApdMap", if APD file has a map type specified and
#   the read map was not given explicitly.
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("readApdUnits", "default", function(filenames, units=NULL, ..., transforms=NULL, cdf=NULL, stratifyBy=c("nothing", "pmmm", "pm", "mm"), addDimnames=FALSE, readMap="byMapType", dropArrayDim=TRUE, verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser");


  apdHeader <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filenames':
  nbrOfArrays <- length(filenames);
  if (nbrOfArrays == 0)
    throw("Cannot read APD units. No APD files specified.");

  filenames <- file.path(dirname(filenames), basename(filenames));
  missing <- filenames[!file.exists(filenames)];
  if (length(missing)) {
    throw("File(s) not found: ", paste(missing, collapse=", "));
  }

  # Argument 'units' and 'cdf':
  if (is.list(cdf) && !is.null(units)) {
    throw("Arguments 'units' must not be specified if argument 'cdf' is a CDF list structure.");
  }

  # Argument 'units':
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- as.integer(units);
    if (any(units < 1))
      stop("Argument 'units' contains non-positive indices.");
  } else {
    stop("Argument 'units' must be numeric or NULL: ", class(units)[1]);
  }

  # Argument 'cdf':
  searchForCdf <- FALSE;
  if (is.null(cdf)) {
    searchForCdf <- TRUE;
  } else if (is.character(cdf)) {
    cdfFile <- file.path(dirname(cdf), basename(cdf));
    if (!file.exists(cdfFile))
      stop("File not found: ", cdfFile);
    cdf <- NULL;
  } else if (is.list(cdf)) {
    aUnit <- cdf[[1]];
    if (!is.list(aUnit))
      throw("Argument 'cdf' is of unknown format: First unit is not a list.");

    groups <- aUnit$groups;
    if (!is.list(groups))
      throw("Argument 'cdf' is of unknown format: Units Does not contain the list 'groups'.");

    # Check for group fields 'indices' or 'x' & 'y' in one of the groups.
    aGroup <- groups[[1]];

    fields <- names(aGroup);
    if ("indices" %in% fields) {
      cdfType <- "indices";
    } else if (all(c("x", "y") %in% fields)) {
      cdfType <- "x";
      searchForCdf <- TRUE;
    } else {
      throw("Argument 'cdf' is of unknown format: The groups contains neither the fields 'indices' nor ('x' and 'y').");
    }
    aUnit <- groups <- aGroup <- NULL; # Not needed anymore
  } else {
    throw("Argument 'cdf' must be a filename, a CDF list structure or NULL: ", mode(cdf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'readMap':
  if (is.null(readMap)) {
    # No probe map specified.
  } else if (identical(readMap, "byMapType")) {
    # Assume that all APD file have the same map type.
    apdHeader <- readApdHeader(filenames[1]);

    mapType <- apdHeader$mapType;
    if (!is.null(mapType)) {
      mapType <- trim(mapType);
      if (is.na(mapType) || nchar(mapType) == 0)
        mapType <- NULL;
    }

    readMap <- mapType;
  }

  if (is.character(readMap)) {
    mapFile <- findApdMap(mapType=readMap);
    if (is.null(mapFile))
      throw("No APD map file found for the given map type: ", readMap);
    readMap <- readApdMap(mapFile)$map;
  }

  # Argument 'dropArrayDim':
  dropArrayDim <- Arguments$getLogical(dropArrayDim);

  # Argument 'addDimnames':
  addDimnames <- Arguments$getLogical(addDimnames);

  # Argument 'transforms':
  if (is.null(transforms)) {
    hasTransforms <- FALSE;
  } else if (is.list(transforms)) {
    if (length(transforms) != nbrOfArrays) {
      throw("Length of argument 'transforms' does not match the number of arrays: ", length(transforms), " != ", nbrOfArrays);
    }
    for (transform in transforms) {
      if (!is.function(transform))
        throw("Argument 'transforms' must be a list of functions.");
    }
    hasTransforms <- TRUE;
  } else {
    throw("Argument 'transforms' must be a list of functions or NULL.");
  }


#  if (length(verbose) != 1)
#    stop("Argument 'verbose' must be a single integer.");
#  verbose <- as.integer(verbose);
#  if (!is.finite(verbose))
#    stop("Argument 'verbose' must be an integer: ", verbose);

  # Argument 'stratifyBy':
  stratifyBy <- match.arg(stratifyBy);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 0. Search for CDF file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (searchForCdf) {
    verbose && enter(verbose, "Searching for CDF file");

    verbose && enter(verbose, "Reading chip type from first APD file (assuming same chip type for all APD files)");
    if (is.null(apdHeader))
      apdHeader <- readApdHeader(filenames[1]);
    chipType <- apdHeader$chipType;
    verbose && exit(verbose);

    verbose && enter(verbose, "Searching for chip type '", chipType, "'");
    cdfFile <- affxparser::findCdf(chipType=chipType);
    if (length(cdfFile) == 0) {
      # If not found, try also where the first APD file is
      opwd <- getwd();
      on.exit(setwd(opwd));
      setwd(dirname(filenames[1]));
      cdfFile <- affxparser::findCdf(chipType=chipType);
      setwd(opwd);
    }
    verbose && exit(verbose);
    if (length(cdfFile) == 0)
      throw("No CDF file for chip type found: ", chipType);

    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 1. Read cell indices for units of interest from the CDF file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(cdf)) {
    verbose && enter(verbose, "Reading cell indices from CDF file");
    cdf <- affxparser::readCdfCellIndices(cdfFile, units=units, stratifyBy=stratifyBy);
    verbose && exit(verbose);

    # Assume 'cdf' contains only "indices" fields.
    indices <- unlist(cdf, use.names=FALSE);
  } else {
    if (cdfType == "indices") {
      # Clean up CDF list structure from other fields than "indices".
      cdfGetFields <- affxparser::cdfGetFields;
      indices <- affxparser::applyCdfGroups(cdf, cdfGetFields, "indices");
      indices <- unlist(cdf, use.names=FALSE);
    } else {
      verbose && enter(verbose, "Calculating cell indices from (x,y) positions");
      verbose && enter(verbose, "Reading chip layout from CDF file");
      cdfHeader <- affxparser::readCdfHeader(cdfFile);
      verbose && exit(verbose);
      cdfGetFields <- affxparser::cdfGetFields;
      x <- unlist(affxparser::applyCdfGroups(cdf, cdfGetFields, "x"), use.names=FALSE);
      y <- unlist(affxparser::applyCdfGroups(cdf, cdfGetFields, "y"), use.names=FALSE);
      ncol <- cdfHeader$cols;
      indices <- as.integer(y * ncol + x + 1);
      x <- y <- ncol <- NULL; # Not needed anymore
      verbose && exit(verbose);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Remapping cell indices?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(readMap)) {
    verbose && enter(verbose, "Remapping cell indices (assuming same map for all files)");
    indices <- readMap[indices];
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Read signals of the cells of interest from the APD file(s)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfCells <- length(indices);
  nbrOfUnits <- length(cdf);

  verbose && enter(verbose, "Reading ", nbrOfUnits, "*", nbrOfCells/nbrOfUnits, "=", nbrOfCells, " cells from ", nbrOfArrays, " APD files");

  for (kk in seq(length=nbrOfArrays)) {
    filename <- filenames[kk];

    verbose && enter(verbose, "Reading APD data for array #", kk);
    apdTmp <- readApd(filename, indices=indices, readMap=NULL, ..., verbose=verbose, .checkArgs=FALSE);
    verbose && exit(verbose);

    # Will future versions of APD hold more than field? /HB 2006-02-28
    dataFields <- setdiff(names(apdTmp), "header");

    if (kk == 1) {
      verbose && enter(verbose, "Allocating return structure");
      # Allocate the return list structure
      apdTmp$header <- NULL;
      apd <- vector("list", length(apdTmp));
      names(apd) <- names(apdTmp);

      for (name in dataFields) {
        apd[[name]] <- matrix(as.double(0), nrow=nbrOfCells, ncol=nbrOfArrays);
      }

      verbose && exit(verbose);
    }

    for (name in dataFields) {
      # Get read signals
      value <- apdTmp[[name]];
      apdTmp[[name]] <- NULL;

      # Transform signals?
      if (hasTransforms && name %in% dataFields) {
        verbose && enter(verbose, "Transform signals for array #", kk);
        value <- transforms[[kk]](value);
        verbose && exit(verbose);
      }
      apd[[name]][,kk] <- value;
    }
  }
  value <- apdTmp <- NULL; # Not needed anymore
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Structure APD data in units and groups according to the CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Structuring data by units and groups");

  # Each group contains the same set of fields
  nbrOfFields <- length(apd);
  fieldNames <- names(apd);
  fields <- vector("list", nbrOfFields);
  names(fields) <- fieldNames;

  # Add a dimension for the arrays, unless only one array is read
  # and the array dimension is not wanted.
  addArrayDim <- (nbrOfArrays >= 2 || !dropArrayDim);

  # For now, assume only one data field.  Extracting the data here, and
  # not in every inner loop, will speed up the whole process about 10 times,
  # e.g. 0.45 seconds instead of 4.3 seconds. /HB 2006-03-04
  data <- apd[[1]];

  seqOfArrays <- list(1:nbrOfArrays);
  offset <- 0;
  res <- lapply(cdf, FUN=function(u) {
    lapply(.subset2(u, "groups"), FUN=function(g) {
      # Same dimensions of all fields
      field <- .subset2(g, 1);  # Faster than g[[1]];
      ncells <- length(field);
      idxs <- offset + 1:ncells;
      offset <<- offset + ncells;

      # Get the target dimension
      dim <- dim(field);
      if (is.null(dim))
        dim <- ncells;

      if (addDimnames) {
        dimnames <- dimnames(field);
        if (is.null(dimnames))
          dimnames <- list(1:ncells);

        # Add an extra dimension for arrays?
        if (addArrayDim) {
          dim <- c(dim, nbrOfArrays);
          dimnames <- c(dimnames, seqOfArrays);
        }

        # Update the data field with dimensions
        # Faster to drop dimensions.
        values <- data[idxs,,drop=TRUE];
        dim(values) <- dim;
        dimnames(values) <- dimnames;
        fields[[1]] <- values;
      } else {
        # Add an extra dimension for arrays?
        if (addArrayDim)
          dim <- c(dim, nbrOfArrays);

        # Update the data field with dimensions
        # Faster to drop dimensions.
        values <- data[idxs,,drop=TRUE];
        if (length(dim) > 1)
          dim(values) <- dim;
        fields[[1]] <- values;
      } # if (addDimnames)

      fields;
    });
  });

  verbose && exit(verbose);

  res;
})


############################################################################
# HISTORY:
# 2006-04-02
# o Using readCdfCellIndices() instead of readCdfUnits().
# 2006-03-28
# o Removed argument 'indexOffset'.
# 2006-03-18
# o Added argument 'indexOffset' and made it one by default (as in R).
# 2006-03-14
# o Added support to obtain the probe map from the map type specified in
#   the APD header.
# 2006-03-04
# o Speed up: By assuming that there is only one single data field in APD
#   files, we can speed up things 10 times.  This is because we do not have
#   to query a list structure in every single inner loop.
# o Removed all gc(). They slow down quite a bit.
# 2006-02-27
# o Created by HB.
############################################################################
