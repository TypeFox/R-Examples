###########################################################################/**
# @RdocClass CnChipEffectFile
#
# @title "The CnChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in a copy-number probe-level
#  models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpChipEffectFile".}
#   \item{combineAlleles}{A @logical indicating if the signals from allele A
#     and allele B are combined or not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically part of a @see "CnChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("CnChipEffectFile", function(..., combineAlleles=FALSE) {
  this <- extend(SnpChipEffectFile(...), c("CnChipEffectFile", uses("CopyNumberDataFile")),
    combineAlleles = combineAlleles
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})

setMethodS3("hasAlleleBFractions", "CnChipEffectFile", function(this, ...) {
  res <- (!this$combineAlleles);
  res;
})

setMethodS3("hasStrandiness", "CnChipEffectFile", function(this, ...) {
  res <- (!this$mergeStrands);
  res;
})


setMethodS3("getParameters", "CnChipEffectFile", function(this, ...) {
  params <- NextMethod("getParameters");
  params$combineAlleles <- this$combineAlleles;
  params;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Retrieves tree list of cell indices for a set of units"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{units}{Unit indices to be retrieved.
#               If @NULL, all units are considered.}
#  \item{...}{Additional arguments passed to \code{getCellIndices()}
#             of @see "ChipEffectFile".}
#  \item{unlist}{If @TRUE, the cell indices are returned as a @vector.}
#  \item{force}{If @TRUE, the cell indices are re-read regardless whether
#     they are already cached in memory or not.}
#  \item{.cache}{(internal) If @TRUE, the result is cached in memory.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @list structure, where each element corresponds to a unit.
#   If argument \code{unlist=TRUE} is passed, an @integer @vector is returned.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getCellIndices", "CnChipEffectFile", function(this, units=NULL, ..., unlist=FALSE, force=FALSE, .cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  # Argument 'units':
  if (is.null(units)) {
  } else if (is.numeric(units)) {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Argument 'unlist':
  unlist <- Arguments$getLogical(unlist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "getCellIndices() for CnChipEffectFile");


  # Supported case?
  combineAlleles <- this$combineAlleles;
  if (unlist && combineAlleles) {
    throw("Unsupported request: Argument 'unlist' have to be TRUE when parameter 'combineAlleles' is TRUE: ", unlist);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force || .cache) {
    chipType <- getChipType(getCdf(this));
    params <- getParameters(this);
    key <- list(method="getCellIndices", class=class(this)[1L],
                chipType=chipType, params=params, units=units, unlist=unlist, ...);
    if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
      key <- getCacheKey(cdf, method="getCellIndices", class=class(this)[1L], chipType=chipType, params=params, units=units, unlist=unlist, ...);
    }
    dirs <- c("aroma.affymetrix", chipType);
    id <- getChecksum(key);
  }

  if (!force) {
    # In memory?
    res <- this$.cellIndices[[id]];
    # On file?
    if (is.null(res)) {
      res <- loadCache(key=list(id), dirs=dirs);
      if (!is.null(res))
        where <- "on file";
    } else {
      where <- "in memory";
    }
    if (!is.null(res)) {
      size <- object.size(res);
      verbose && printf(verbose, "Returning value cached %s: %.1fMB\n",
                                                   where, size/1024^2);
      verbose && exit(verbose);
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get units in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(cdf));
  }
  nbrOfUnits <- length(units);

  cells <- lapplyInChunks(units, function(unitChunk) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Get and restructure cell indices
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## NOTE: NextMethod() does not work from within another function
##    cells <- NextMethod("getCellIndices", units=unitChunk, force=force, .cache=FALSE, verbose=verbose);
    cells <- getCellIndices.SnpChipEffectFile(this, units=unitChunk, ...,
              unlist=unlist, force=force, .cache=FALSE, verbose=verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Combine alleles?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # If combining alleles, return only every second group.
    # In order to improve readability we merge the names of alleles groups
    # combined, e.g. groups 'C' and 'G' become group 'CG'.
    if (combineAlleles) {
      verbose && enter(verbose, "Combining alleles");
      # Hard-wiring 1, 2 & 4 groups speed things up 3 times!

      cells <- .applyCdfGroups(cells, function(groups) {
        ngroups <- length(groups);
        names <- names(groups);
        if (ngroups == 4L) {
          groups <- .subset(groups, c(1L,3L));
          names <- paste(.subset(names, c(1L,3L)), .subset(names, c(2L,4L)), sep="");
        } else if (ngroups == 2L) {
          groups <- .subset(groups, 1L);
          names <- paste(.subset(names, 1L), .subset(names, 2L), sep="");
        } else if (ngroups == 1L) {
          groups <- .subset(groups, 1L);
        } else {
          odds <- seq(from=1L, to=ngroups, by=2L);
          evens <- seq(from=2L, to=ngroups, by=2L);
          groups <- .subset(groups, odds);
          names <- paste(.subset(names, odds), .subset(names, evens), sep="");
        }
        names(groups) <- names;
        groups;
      }) # applyCdfGroups()
      verbose && printf(verbose, "Number of units: %d\n", length(cells));
      verbose && exit(verbose);
    } # if (combineAlleles)

    cells;
  }, chunkSize=100e3, verbose=less(verbose)) # lapplyInChunks()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Unlist?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (unlist) {
    cells <- unlist(cells, use.names=FALSE);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (.cache) {
    # In-memory or on-file cache?
    size <- object.size(cells);
      verbose && printf(verbose, "Object size: %.1fMB\n", size/1024^2);
    if (size < 10e6) {
      # In-memory cache for objects < 10Mb.
      this$.cellIndices <- list();
      this$.cellIndices[[id]] <- cells;
      verbose && cat(verbose, "Result cached in memory");
    } else {
      # On-file cache
      # Keep, in-memory cache.
      if (!is.list(this$.cellIndices))
        this$.cellIndices <- list();
      this$.cellIndices[[id]] <- NULL;
      saveCache(cells, key=list(id), dirs=dirs);
      verbose && cat(verbose, "Result cached to file");
    }
  }

  verbose && exit(verbose);

  cells;
}, protected=TRUE) # getCellIndices()


setMethodS3("readUnits", "CnChipEffectFile", function(this, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Check for cached data
  key <- list(method="readUnits", class=class(this)[1],
              pathname=getPathname(this),
              combineAlleles=this$combineAlleles, ...);
  id <- getChecksum(key);
  res <- this$.readUnitsCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.CnChipEffectFile(): Returning cached data");
    return(res);
  }

  # Retrieve the data
  res <- NextMethod("readUnits", force=TRUE, cache=FALSE, verbose=verbose);


  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.CnChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
  }

  res;
})


setMethodS3("mergeStrands", "CnChipEffectFile", function(this, ...) {
  cfM <- NextMethod("mergeStrands");
  cfM$combineAlleles <- this$combineAlleles;
  cfM;
})


setMethodS3("getNumberOfFilesAveraged", "CnChipEffectFile", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Querying ", class(this)[1], " for the number of arrays used to calculate average");
  pathname <- getPathname(this);
  verbose && cat(verbose, "Pathname: ", pathname);

  # Make sure it is a file calculated by getAverageFile()
  if (!isAverageFile(this)) {
    throw("Cannot retrieve the number of arrays used when averaging. The file is not generated by getAverageFile(), because the filename does not start with '.average-': ", pathname);
  }

  # Get (unit, group, cell) map
  ugcMap <- getUnitGroupCellMap(this, ..., verbose=less(verbose, 5));

  # Keep only first group (in case mergeStrands or combineAlleles is FALSE)
  # To please R CMD check
  group <- NULL; rm(list="group");
  ugcMap <- subset(ugcMap, group == 1);

  cells <- ugcMap[,"cell", drop=TRUE];
  # Not needed anymore
  ugcMap <- NULL;

  verbose && cat(verbose, "Cell indices:");
  verbose && str(verbose, cells);

  verbose && enter(verbose, "Reading data");
  data <- .readCel(pathname, indices=cells, readIntensities=FALSE,
                  readPixels=TRUE);
  verbose && exit(verbose);

  ns <- data$pixels;
  # Not needed anymore
  data <- NULL;

  verbose && exit(verbose);

  ns;
}, protected=TRUE); # getNumberOfFilesAveraged()


############################################################################
# HISTORY:
# 2012-11-29
# o ROBUSTNESS: Added protection for getCellIndices(..., unlist=TRUE)
#   for SnpChipEffectFile and CnChipEffectFile.
# o Added Rdoc help for getCellIndices() for CnChipEffectFile.
# 2009-11-18
# o Added getNumberOfFilesAveraged() for CnChipEffectFile.  It should only
#   be used for files generated by getAverageFile(), otherwise an error
#   is thrown.
# 2009-11-17
# o Now CnChipEffectFile implements CopyNumberDataFile.
# 2008-02-19
# o BUG FIX: getCellIndices() of CnChipEffectFile would return an error
#   if 'units==NULL'.
# 2007-03-01
# o BUG FIX: getCellIndices() would give "Error in fcn(.subset2(unit,
#   "groups"), ...) : object "odds" not found" for units with other than
#   1, 2, or 4 groups.
# 2007-01-20
# o Added mergeStrands().
# 2007-01-06
# o Now getCellIndices() caches large objects to file and small in memory.
# o Made getCellIndices() three times faster by some hardwired code.
# 2006-11-28
# o Added readUnits() to override caching mechanism of superclasses.
# 2006-09-20
# o BUG FIX: Typo. Remove an argument but tried to use inside.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-12
# o Updated and probably working. When combining alleles, the names of the
#   groups returned consist of the allele A and allele group names.
# 2006-09-11
# o Created.
############################################################################
