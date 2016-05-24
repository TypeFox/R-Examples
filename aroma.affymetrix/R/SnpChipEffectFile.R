###########################################################################/**
# @RdocClass SnpChipEffectFile
#
# @title "The SnpChipEffectFile class"
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
#   \item{mergeStrands}{Specifies if the strands are merged or not for these
#      estimates.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically part of a @see "SnpChipEffectSet".
# }
#*/###########################################################################
setConstructorS3("SnpChipEffectFile", function(..., mergeStrands=FALSE) {
  this <- extend(ChipEffectFile(...), "SnpChipEffectFile",
    "cached:.cellIndices" = NULL,
    mergeStrands = mergeStrands
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getParameters", "SnpChipEffectFile", function(this, ...) {
  params <- NextMethod("getParameters");
  params$mergeStrands <- this$mergeStrands;
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
setMethodS3("getCellIndices", "SnpChipEffectFile", function(this, units=NULL, ..., unlist=FALSE, force=FALSE, .cache=TRUE, verbose=FALSE) {
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

  verbose && enter(verbose, "getCellIndices() for SnpChipEffectFile");

  # Supported case?
  mergeStrands <- this$mergeStrands;
  if (unlist && mergeStrands) {
    throw("Unsupported request: Argument 'unlist' have to be TRUE when parameter 'mergeStrands' is TRUE: ", unlist);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!force || .cache) {
    chipType <- getChipType(cdf);
    params <- getParameters(this);
    key <- list(method="getCellIndices", class=class(this)[1L],
                pathname=getPathname(this),  ## <= WRONG! /HB 2012-11-29
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
    cells <- getCellIndices.ChipEffectFile(this, units=unitChunk, ...,
              unlist=unlist, force=force, .cache=FALSE, verbose=verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Merge strands?
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # If merging strands, we only need half the number of chip-effect
    # parameters per unit group.  Example:
    # (a) mergeStrands=FALSE:
    #    Fit by strand and allele:        #groups=4, #chip effects=4
    #    (same but single-stranded SNP)   #groups=2, #chip effects=2
    # (b) mergeStrands=TRUE:
    #    Merge strands, fit by allele:    #groups=4, #chip effects=2
    #    (same but single-stranded SNP)   #groups=2, #chip effects=2
    if (mergeStrands) {
      verbose && enter(verbose, "Merging strands");
      cells <- .applyCdfGroups(cells, function(groups) {
        ngroups <- length(groups);
        if (ngroups == 4L) {
          .subset(groups, c(1L,2L));
        } else if (ngroups == 2L) {
          .subset(groups, c(1L,2L));
        } else if (ngroups == 1L) {
          .subset(groups, 1L);
        } else if (ngroups > 4L && ngroups %% 2L == 0L) {
          .subset(groups, c(1L,2L));
        } else {
          # groups[1:ceiling(ngroups/2L)];
          # groups[1:round((ngroups+1L)/2L)];
          .subset(groups, 1:round((ngroups+1L)/2L));
        }
      }) # applyCdfGroups()
      verbose && printf(verbose, "Number of units: %d\n", length(cells));
      verbose && exit(verbose);
    } # if (this$mergeStrands)

    cells;
  }, chunkSize=100e3, verbose=less(verbose));


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


setMethodS3("mergeStrands", "SnpChipEffectFile", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mergeStrandsMatrix <- function(y, ...) {
    n <- nrow(y);
    if (n == 2) {
      y[1,] <- colMeans(y[1:2,,drop=FALSE], na.rm=TRUE);
      y[2,] <- NA;
      attr(y, "groups") <- 1;
    } else if (n == 4) {
      y[1,] <- colMeans(y[1:2,,drop=FALSE], na.rm=TRUE);
      y[3,] <- colMeans(y[3:4,,drop=FALSE], na.rm=TRUE);
      y[c(2,4),] <- NA;
    } else if (n > 4 && n %% 2 == 0) {
      # All other even-numbered groups
      odd <- seq(from=1, to=n, by=2);
      for (idx in odd) {
        y[idx,] <- colMeans(y[idx:(idx+1),,drop=FALSE], na.rm=TRUE);
      }
      even <- seq(from=2, to=n, by=2);
      y[even,] <- NA;
    }
    y;
  }

  if (this$mergeStrands) {
    throw("Strands have already been merged for this ", class(this)[1]);
  }

  cfM <- mergeGroups(this, fcn=mergeStrandsMatrix, ...);
  cfM$mergeStrands <- TRUE;

  cfM;
}, protected=TRUE)



setMethodS3("readUnits", "SnpChipEffectFile", function(this, ..., force=FALSE, cache=FALSE, verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Check for cached data
  key <- list(method="readUnits", class=class(this)[1],
              mergeStrands=this$mergeStrands, ...);
  id <- getChecksum(key);
  res <- this$.readUnitsCache[[id]];
  if (!force && !is.null(res)) {
    verbose && cat(verbose, "readUnits.SnpChipEffectFile(): Returning cached data");
    return(res);
  }

  # Retrieve the data
  res <- NextMethod("readUnits", force=TRUE, cache=FALSE, verbose=verbose);

  # Store read units in cache?
  if (cache) {
    verbose && cat(verbose, "readUnits.SnpChipEffectFile(): Updating cache");
    this$.readUnitsCache <- list();
    this$.readUnitsCache[[id]] <- res;
  }

  res;
})


############################################################################
# HISTORY:
# 2012-11-29
# o SPEEDUP: Improved the caching mechanism for getCellIndices() for
#   SnpChipEffectFile and CnChipEffectFile.
# o ROBUSTNESS: Added protection for getCellIndices(..., unlist=TRUE)
#   for SnpChipEffectFile and CnChipEffectFile.
# o Added Rdoc help for getCellIndices() for SnpChipEffectFile.
# 2008-02-22
# o Updated getCellIndices() and mergeStrands() to handle any even-numbered
#   unit groups beyond two and four groups.
# 2007-01-20
# o Added mergeStrands().
# 2007-01-07
# o Now getCellIndices() caches large objects to file and small in memory.
# 2006-12-18
# o BUG FIX: getCellIndices() would return a single group instead of two,
#   for a single-stranded SNP when mergeStrands=TRUE.  See for instance
#   unit SNP_A-1780520 in the Mapping250K_Nsp chip.
# 2006-11-28
# o Added readUnits() to override caching mechanism of superclasses.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-12
# o Updated.  Now the names of the groups reflects the allele names as
#   expected.
# 2006-09-11
# o Created.
############################################################################
