setMethodS3("nbrOfCells", "AffymetrixCdfFile", function(this, ...) {
  as.integer(prod(getDimension(this, ...)));
})

setMethodS3("nbrOfUnits", "AffymetrixCdfFile", function(this, ...) {
  getHeader(this)$probesets;
})

setMethodS3("nbrOfQcUnits", "AffymetrixCdfFile", function(this, ...) {
  getHeader(this)$qcprobesets;
})


###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod nbrOfGroupsPerUnit
#
# @title "Gets the number of groups in each unit"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @vector of @integers.
# }
#
# \details{
#   Once read from file, this information is cached in memory for efficiency.
#   The cache can be cleared by calling \code{gc(cdf)}.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("nbrOfGroupsPerUnit", "AffymetrixCdfFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  sizes <- this$.nbrOfGroupsPerUnit;

  if (force || is.null(sizes)) {
    # Check in file cache
    chipType <- getChipType(this);
    key <- list(method="nbrOfGroupsPerUnit", class=class(this)[1],
                chipType=chipType);
    if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
      key <- getCacheKey(this, method="nbrOfGroupsPerUnit", chipType=chipType);
    }
    dirs <- c("aroma.affymetrix", chipType);
    if (force) {
      sizes <- NULL;
    } else {
      sizes <- loadCache(key=key, dirs=dirs);
    }

    if (is.null(sizes)) {
      verbose && enter(verbose, "Reading number of groups for *all* units");
      sizes <- .readCdfGroupNames(getPathname(this));
      sizes <- restruct(this, sizes, verbose=less(verbose, 5));
      sizes <- lapply(sizes, FUN=length);
      sizes <- unlist(sizes, use.names=FALSE);
      saveCache(sizes, key=key, dirs=dirs);
      verbose && exit(verbose);
    }

    this$.nbrOfGroupsPerUnit <- sizes;
  }

  if (!is.null(units))
    sizes <- sizes[units];

  sizes;
}, private=TRUE)


setMethodS3("nbrOfCellsPerUnitGroup", "AffymetrixCdfFile", function(this, units=NULL, ..., useNames=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  chipType <- getChipType(this);
  verbose && enter(verbose, "Getting the number of cells per unit group");
  verbose && cat(verbose, "Chip type: ", chipType);

  key <- list(method="nbrOfCellsPerUnitGroup", class=class(this)[1],
              chipType=chipType, useNames=useNames);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="nbrOfCellsPerUnitGroup", chipType=chipType, useNames=useNames);
  }
  dirs <- c("aroma.affymetrix", chipType);
  if (force) {
    counts <- NULL;
  } else {
    verbose && enter(verbose, "Checking for cached results");
    counts <- loadCache(key=key, dirs=dirs);
    if (!is.null(counts))
      verbose && cat(verbose, "Cached results found.");
    verbose && exit(verbose);
  }

  if (is.null(counts)) {
    verbose && enter(verbose, "Reading cell counts from CDF");
    counts <- .readCdfNbrOfCellsPerUnitGroup(getPathname(this));
    verbose && exit(verbose);

    counts <- restruct(this, counts, verbose=less(verbose, 5));

    if (!useNames) {
      names(counts) <- NULL;
      counts <- lapply(counts, FUN=function(x) { names(x) <- NULL; x });
    }

    saveCache(counts, key=key, dirs=dirs);
  }

  # Subset?
  if (!is.null(units)) {
    counts <- counts[units];
  }

  verbose && exit(verbose);

  counts;
}, private=TRUE)



setMethodS3("nbrOfCellsPerUnit", "AffymetrixCdfFile", function(this, units=NULL, ..., useNames=FALSE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  chipType <- getChipType(this);
  verbose && enter(verbose, "Getting the number of cells per unit");
  verbose && cat(verbose, "Chip type: ", chipType);

  key <- list(method="nbrOfCellsPerUnit", class=class(this)[1],
              chipType=chipType, useNames=useNames);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="nbrOfCellsPerUnit", chipType=chipType, useNames=useNames);
  }
  dirs <- c("aroma.affymetrix", chipType);
  if (force) {
    counts <- NULL;
  } else {
    verbose && enter(verbose, "Checking for cached results");
    counts <- loadCache(key=key, dirs=dirs);
    if (!is.null(counts))
      verbose && cat(verbose, "Cached results found.");
    verbose && exit(verbose);
  }

  if (is.null(counts)) {
    verbose && enter(verbose, "Getting number of cells per unit group");
    counts <- nbrOfCellsPerUnitGroup(this, units=NULL, useNames=useNames,
                                           ..., verbose=less(verbose, 1));

    verbose && enter(verbose, "Summing per unit");
    # Sum per unit
    counts <- sapply(counts, FUN=sum);
    verbose && exit(verbose);

    verbose && exit(verbose);

    saveCache(counts, key=key, dirs=dirs);
  }

  # Subset?
  if (!is.null(units)) {
    counts <- counts[units];
  }

  verbose && exit(verbose);

  counts;
}, private=TRUE)



############################################################################
# HISTORY:
# 2008-10-09
# o Added nbrOfCellsPerUnit() and nbrOfCellsPerUnitGroup().
# o Added verbose output to internal restruct().  Is that ever used?!?
# o Renamed getUnitSizes() to nbrOfGroupsPerUnit().
# 2008-09-06
# o BUG FIX: getUnitTypes() of AffymetrixCdfFile did not return a
#   name map for the unit types if a subset was units was selected.
# 2008-08-09
# o BUG FIX: getUnitTypes() of AffymetrixCdfFile would not return the
#   correct integer for binary CDFs.  Now it uses readCdf() instead of
#   readCdfUnits() of affxparser.
# 2008-07-26
# o Now the integer vector returned by getUnitTypes() also has an attribute
#   'types' explain what the different values are.
# 2008-07-23
# o Now getGenomeInformation() and getSnpInformation() reports the reason
#   for why it thinks the located object is incompatible with the CDF.
# 2008-05-18
# o Now AffymetrixCdfFile "provides" the UnitNamesInterface.
# 2008-05-09
# o Now inherits from AromaChipTypeAnnotationFile.
# 2008-04-12
# o BUG FIX: getChipType(..., fullname=FALSE) would return the chip type as
#   the 'tags' attribute if there were no tags.
# 2008-04-08
# o Added byName().
# 2008-03-11
# o BUG FIX: Calling readUnits() of an AffymetrixCdfFile without specifying
#   the 'units' argument gave an error.  Thanks Tim Keighley, CSIRO, Sydney
#   for reporting this.
# 2008-02-21
# o Added getGroupDirections() for AffymetrixCdfFile.
# o Added getUnitTypes() for AffymetrixCdfFile.
# 2008-01-20
# o Now getSnpInformation() searches for UFL files with tags.
# 2008-01-19
# o Now getGenomeInformation() searches for UGP files with tags.
# 2007-12-09
# o Now get- and set- Genome/SnpInformation() asserts that the annotation
#   file objects are compatible with the CDF.  At least for UGP & UFL files.
# 2007-12-08
# o Added setGenomeInformation() & setSnpInformation() to AffymetrixCdfFile.
# o Now construct AffymetrixCdfFile$fromName("HuEx-1_0-st-v2", tags="core")
#   can be used to locate 'HuEx-1_0-st-v2,core.CDF'.
# 2007-09-10
# o Now getGenomeInformation() of AffymetrixCdfFile recognizes UGP files
#   as well and before dChip genome information files.
# 2007-09-06
# o Now identifyCells() utilized the below to save memory.
# o Added argument 'useNames=TRUE' and 'unlist=FALSE' to getCellIndicies()
#   of AffymetrixCdfFile.  These can be use to save memory.
# 2007-08-17
# o Made getCellIndices() of AffymetrixCdfFile more memory effiencent by
#   reading and transforming data in chunks.
# 2007-08-09
# o Now convertCdf() generates a CDF file upper-case extension *.CDF.
# 2007-08-02
# o Renamed fromChipType() of AffymetrixCdfFile to byChipType().
# 2007-07-09
# o Added getFileFormat() to AffymetrixCdfFile.  This is also reported
#   by the print() method.
# 2007-03-28
# o Added argument 'cache=TRUE' to getCellIndices().
# 2007-03-26
# o Added a few more gc().
# o BUG FIX: isPm() did not work when querying a subset of the units.
# 2007-02-22
# o Now findByChipType() recognizes Windows shortcuts too.
# 2007-02-21
# o Now findByChipType() passes '...' to underlying function.
# 2007-02-14
# o BUG FIX: When "tagifying" monocell, getSnpInformation() and
#   getGenomeInformation() was looking for the incorrect chip type.
# 2007-02-12
# o Added argument 'main' to getChipType().
# 2007-02-08
# o Now findByChipType() handles monocell CDFs specially; monocell CDFs can
#   still be put in the same directory as the parent CDF.
# 2007-02-06
# o Added findByChipType().
# 2007-01-16
# o Now all cache keys contains method name, class name, and chip type.
# 2007-01-10
# o Reordered internally in createMonoCell() preparing for code to read
#   *and* write monocell CDFs in chunks.  It should not be too hard.
#   We need to update affxparser with writeCdfHeader(), writeCdfQcUnits()
#   and writeCdfUnits(), and are basically already in there, but as
#   private functions.
# o Removed some unnecessary group fields in 'destUnits' saving us approx
#   220-230Mb out of 1.9GB for the Mapping250K_Nsp.  Still need to find
#   another solution to get down below 1GB. One thing that takes up a lot
#   of memory is that the unit and group directions are stored as character
#   strings and not integers.
# o BUG FIX: The most recent createMonoCell() would create CDFs with all
#   cell indices being equal to one.  Added more verbose output and some
#   garbage collection to this function too.
# 2007-01-06
# o Added argument 'force' to getCellIndices().
# o Now getCellIndices() only caches object < 10 MB RAM.
# o Optimized identifyCells() to only cache data in rare cases.
# 2006-12-18 /KS
# o Made global replacement "block" -> "group".
# 2006-12-14
# o Added convertUnits().
# 2006-09-27
# o Now fromFile() tries to create an instance of the subclasses (bottom up)
#   first.  This will make it possible to automatically define SNP CDFs.
# 2006-09-26
# o Now getGenomeInformation() and getSnpInformation() ignores suffices of
#   the chip-type string. This makes it possible to retrive annotation data
#   also for custom chips.
# 2006-09-17
# o Added an in-memory cache for getCellIndices().
# 2006-09-16
# o Added getGenomeInformation() and stextChipType().
# 2006-09-14
# o BUG FIX: Fractional value of 'indices' of identifyCells().
# 2006-09-10
# o BUG FIX: createMonoCell() where resetting the cell indices for each
#   chunk.
# o Simple benchmarking of createMonoCell(): IBM Thinkpad A31 1.8GHz 1GB:
#   Mapping50K_Hind to mono cells CDF takes ~13 mins.  Again, it is
#   writeCdf() that is slow.  KH is working on improving this.
# 2006-09-08
# o Added equals() to compare to CDF object.
# o Added convert() to convert a CDF into another version by convertCdf().
# 2006-08-25
# o Added getFirstCellIndices().  This may be used by methods to identify
#   unit or unit groups for which no probe signals have been assigned yet.
# 2006-08-24
# o Added reconstruct().
# o Added Rdoc comments.
# 2006-08-19
# o Added getUnitSizes().  Useful in order to store parameter estimates
#   of probeset-summary models.
# 2006-08-11
# o Created.
############################################################################
