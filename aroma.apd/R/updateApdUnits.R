#########################################################################/**
# @RdocDefault updateApdUnits
#
# @title "Updates an Affymetrix probe data (APD) file by units (probesets)"
#
# @synopsis
#
# \description{
#   @get "title" by using the unit and group definitions in the
#   corresponding Affymetrix CDF file.
# }
#
# \arguments{
#   \item{filename}{The filename of the APD file.}
#   \item{units}{An @integer @vector of unit indices specifying which
#     units to be read.  If @NULL, all units are updated.}
#   \item{data}{A @numeric @vector of data elements to be assigned.}
#   \item{...}{Additional arguments passed to @see "updateApd",
#     e.g. \code{writeMap}.}
#   \item{cdf}{A @character filename of a CDF file, or a CDF @list
#     structure.  If @NULL, the CDF file is searched for by
#     @see "affxparser::findCdf".}
#   \item{stratifyBy}{Argument passed to low-level method
#     @see "affxparser::readCdfCellIndices".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \examples{\dontrun{#See ?createApd for an example.}}
#
# \seealso{
#   @see "readApdUnits" to read unit by units.
#   @see "updateApd" to update cell by cell.
# }
#
# @keyword "file"
# @keyword "IO"
#*/#########################################################################
setMethodS3("updateApdUnits", "default", function(filename, units=NULL, data, ..., cdf=NULL, stratifyBy=c("nothing", "pmmm", "pm", "mm"), verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  filename <- Arguments$getReadablePathname(filename, mustExist=TRUE);

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

  # Argument 'data':
  if (!is.list(data)) {
    throw("Argument 'data' is not a list: ", mode(data));
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


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 0. Search for CDF file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (searchForCdf) {
    verbose && enter(verbose, "Searching for CDF file");

    verbose && enter(verbose, "Reading chip type from first APD file");
    apdHeader <- readApdHeader(filename);
    chipType <- apdHeader$chipType;
    verbose && exit(verbose);

    verbose && enter(verbose, "Searching for chip type '", chipType, "'");
    cdfFile <- affxparser::findCdf(chipType=chipType);
    if (length(cdfFile) == 0) {
      # If not found, try also where the first APD file is
      opwd <- getwd();
      on.exit(setwd(opwd));
      setwd(dirname(filename));
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
    cdfType <- "flat";
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 2. Assert that 'data' has the same structure as 'cdf'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nbrOfUnits <- length(cdf);
  if (nbrOfUnits != length(data)) {
    throw("Argument 'data' and the CDF structure are of different lengths: ",
                                           nbrOfUnits, " != ", length(data));
  }

  # For each CDF unit
  for (uu in seq(length=nbrOfUnits)) {
    unit <- .subset2(cdf, uu);  # Faster than cdf[[uu]]
    cdfGroups <- .subset2(unit, "groups");
    dataUnit <- .subset2(data, uu);
    nbrOfGroups <- length(cdfGroups);
    if (nbrOfGroups != length(dataUnit)) {
      unitName <- names(cdf)[uu];
      throw("The number of groups in unit ", unitName, " does not match the corresponding data unit: ", length(cdfGroups), " != ", length(dataUnit));
    }

    if (nbrOfGroups == 0)
      next;

    # For unit group
    for (gg in 1:nbrOfGroups) {
      cdfField <- .subset2(.subset2(cdfGroups, gg), 1);
      dataField <- .subset2(.subset2(dataUnit, gg), 1);
      if (length(cdfField) != length(dataField)) {
        unitName <- names(cdf)[uu];
        throw("The number of values in field #", gg, " in the groups of unit ", unitName, " does not match the corresponding data unit field: ", length(cdfField), " != ", length(dataField));
      }
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 3. Extract cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cdfType == "flat") {
    indices <- unlist(cdf, use.names=FALSE);
  } else if (cdfType == "indices") {
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
    x <- unlist(affxparser::applyCdfGroups(cdf, cdfGetFields, "x"), use.names=TRUE);
    y <- unlist(affxparser::applyCdfGroups(cdf, cdfGetFields, "y"), use.names=TRUE);
    ncol <- cdfHeader$cols;
    indices <- as.integer(y * ncol + x + 1);
    x <- y <- ncol <- NULL; # Not needed anymore
    verbose && exit(verbose);
  }

  # Flatten 'data'
  data <- unlist(data, use.names=FALSE);

  # Assert same length again
  if (length(indices) != length(data)) {
    throw("The number of cells and the number of elements in data does not match: ", length(indices), " != ", length(data));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  updateApd(filename, indices=indices, data=data, ..., .checkArgs=FALSE);
})


############################################################################
# HISTORY:
# 2006-06-18
# o Removed argument 'addDimnames', which was never used.
# 2006-04-02
# o Using readCdfCellIndices() instead of readCdfUnits().
# 2006-03-28
# o Removed argument 'indexOffset' and 'writeMap'.
# 2006-03-18
# o Added argument 'indexOffset' and made it one by default (as in R).
# 2006-03-17
# o Created.
############################################################################
