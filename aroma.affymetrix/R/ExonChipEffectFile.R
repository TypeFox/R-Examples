###########################################################################/**
# @RdocClass ExonChipEffectFile
#
# @title "The ExonChipEffectFile class"
#
# \description{
#  @classhierarchy
#
#  This class represents estimates of chip effects in the probe-level models.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ChipEffectFile".}
#   \item{mergeGroups}{Specifies if the groups are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "KS, HB"
#
# \seealso{
#   An object of this class is typically part of a @see "ExonChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("ExonChipEffectFile", function(..., mergeGroups=FALSE) {
  this <- extend(ChipEffectFile(...), "ExonChipEffectFile",
    "cached:.cellIndices" = NULL,
    mergeGroups = mergeGroups
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getParameters", "ExonChipEffectFile", function(this, ...) {
  params <- NextMethod("getParameters");
  params$mergeGroups <- this$mergeGroups;
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
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("getCellIndices", "ExonChipEffectFile", function(this, ..., unlist=FALSE, force=FALSE, .cache=TRUE, verbose=FALSE) {
  # Argument 'unlist':
  unlist <- Arguments$getLogical(unlist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "getCellIndices() for ExonChipEffectFile");

  # Supported case?
  mergeGroups <- this$mergeGroups;
  if (unlist && mergeGroups) {
    throw("Unsupported request: Argument 'unlist' have to be TRUE when parameter 'mergeGroups' is TRUE: ", unlist);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force || .cache) {
    cdf <- getCdf(this);
    chipType <- getChipType(cdf);
    params <- getParameters(this);
    key <- list(method="getCellIndices", class=class(this)[1],
                pathname=getPathname(this),
                chipType=chipType, params=params, unlist=unlist, ...);
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
  # Get and restructure cell indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cells <- NextMethod("getCellIndices", .cache=FALSE, verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Merge groups?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # If merging groups, we only need one chip-effect parameter per unit
  if (mergeGroups) {
    verbose && enter(verbose, "Merging groups");
    cells <- .applyCdfGroups(cells, function(groups) {
      .subset(groups, 1L);
    })
    verbose && exit(verbose);
  }


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
    if (object.size(cells) < 10e6) {
      # In-memory cache for objects < 10Mb.
      this$.cellIndices <- list();
      this$.cellIndices[[id]] <- cells;
      verbose && cat(verbose, "Cached in memory");
    } else {
      # On-file cache
      # Keep, in-memory cache.
      if (!is.list(this$.cellIndices))
        this$.cellIndices <- list();
      this$.cellIndices[[id]] <- NULL;
      saveCache(cells, key=list(id), dirs=dirs);
      verbose && cat(verbose, "Cached to file");
    }
  }

  verbose && exit(verbose);

  cells;
}, protected=TRUE) # getCellIndices()


setMethodS3("readUnits", "ExonChipEffectFile", function(this, ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Check for cached data
  key <- list(method="readUnits", class=class(this)[1],
              mergeGroups=this$mergeGroups, ...);
  id <- getChecksum(key);
  res <- this$.readUnitsCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.ExonChipEffectFile(): Returning cached data");
    return(res);
  }

  # Retrieve the data
  res <- NextMethod("readUnits", force=TRUE, cache=FALSE, verbose=verbose);

  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.ExonChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
  }

  res;
})


## setMethodS3("findUnitsTodo", "ExonChipEffectFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
##   # Argument 'verbose':
##   verbose <- Arguments$getVerbose(verbose);
##   if (verbose) {
##     pushState(verbose);
##     on.exit(popState(verbose));
##   }
##
##
##   verbose && enter(verbose, "Identifying non-fitted units in chip-effect file");
##   verbose && cat(verbose, "Pathname: ", getPathname(this));
##
##
##   idxs <- NULL;
##   if (is.null(units)) {
##     # Look up chip-type and parameter specific but data set independent data
##     cdf <- getCdf(this);
##     chipType <- getChipType(cdf);
##     key <- list(method="findUnitsTodo", class=class(this)[1],
##                 chipType=chipType, params=getParameters(this));
##     dirs <- c("aroma.affymetrix", chipType);
##     if (!force) {
##       idxs <- loadCache(key, dirs=dirs);
##       if (!is.null(idxs))
##         verbose && cat(verbose, "Found indices cached on file");
##     }
##   }
##
##   if (is.null(idxs)) {
##     verbose && enter(verbose, "Identifying CDF units");
##
##     verbose && enter(verbose, "Reading CDF cell indices");
##     idxs <- getCellIndices(this, units=units, verbose=less(verbose));
##     # Example:
##     #  $ 2315554:List of 1
##     #   ..$ groups:List of 1
##     #   .. ..$ groupA:List of 1
##     #   .. .. ..$ indices: int 1
##     #   .. ..$ groupB:List of 1
##     #   .. .. ..$ indices: int 23
##     verbose && exit(verbose);
##
##     verbose && enter(verbose, "Extracting first CDF block for each unit");
##     idxs <- applyCdfGroups(idxs, .subset2, 1);
##     # Example:
##     #  $ 2315554:List of 1
##     #   ..$ groups:List of 1
##     #   .. ..$ indices: int 1
##     ## The below makes no difference.  Thus, this function give exactly
##     ## the same result as the one in the super class. Note, mergeGroups
##     ## has already been taken care of inside getCellIndices() /HB 2007-08-17
##     if (this$mergeGroups) {
##       # do a further round of reduction
##       idxs <- applyCdfGroups(idxs, .subset2, 1);
##       # Example:
##       #  $ 2315554:List of 1
##       #   ..$ groups: int 1
##     }
##     verbose && exit(verbose);
##
##     idxs <- unlist(idxs, use.names=FALSE);
##
##     if (is.null(units)) {
##       verbose && enter(verbose, "Saving to file cache");
##       saveCache(idxs, key=key, dirs=dirs);
##       verbose && exit(verbose);
##     }
##
##     verbose && exit(verbose);
##   }
##
##
##   # Read one cell from each unit
##   verbose && enter(verbose, "Reading data for these ", length(idxs), " cells");
##   value <- readCel(getPathname(this), indices=idxs, readIntensities=FALSE,
##                    readStdvs=TRUE, readPixels=FALSE)$stdvs;
##   verbose && exit(verbose);
##
##
##   # Identify units for which the stdvs <= 0.
##   value <- which(value <= 0);
##
##   if (!is.null(units))
##     value <- units[value];
##   verbose && cat(verbose, "Looking for stdvs <= 0 indicating non-estimated units:");
##   verbose && str(verbose, value);
##
##   verbose && exit(verbose);
##
##   value;
## })


############################################################################
# HISTORY:
# 2007-08-17 /HB
# o Removed findUnitsTodo() from ExonChipEffectFile, because it gave the
#   same result as in the one in superclass ChipEffectFile.
# 2007-02-08 /KS
# o Created (based on SnpChipEffectFile.R following chat with HB on
#   2007-02-07).
############################################################################
