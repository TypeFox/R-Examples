###########################################################################/**
# @RdocClass AffymetrixCdfFile
#
# @title "The AffymetrixCdfFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixCdfFile object represents a generic Affymetrix CDF file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaChipTypeAnnotationFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB, KS"
#*/###########################################################################
setConstructorS3("AffymetrixCdfFile", function(...) {
  this <- extend(AromaChipTypeAnnotationFile(...), c("AffymetrixCdfFile",
            uses("UnitNamesFile", "UnitTypesFile", "AromaPlatformInterface")),
    "cached:.header" = NULL,
    "cached:.unitNames" = NULL,
    "cached:.unitTypes" = NULL,
    "cached:.cellIndices" = NULL,
    "cached:.isPm" = NULL,
    "cached:.gi" = NULL,
    "cached:.si" = NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getExtensionPattern", "AffymetrixCdfFile", function(static, ...) {
  "[.](cdf|CDF)$";
}, static=TRUE, protected=TRUE)



setMethodS3("getUnitNamesFile", "AffymetrixCdfFile", function(this, ...) {
  this;
}, protected=TRUE)


setMethodS3("getUnitTypesFile", "AffymetrixCdfFile", function(this, ...) {
  this;
}, protected=TRUE)


setMethodS3("getFileFormat", "AffymetrixCdfFile", function(this, ...) {
  pathname <- getPathname(this);
  if (!isFile(pathname)) return(NA_character_)

  # Read CDF header
  raw <- readBin(pathname, what=raw(), n=10);

  if (raw[1] == 59)
    return("v5 (binary; CC)");

  if (raw[1] == 67)
    return("v4 (binary; XDA)");

  if (rawToChar(raw[1:5]) == "[CDF]")
    return("v3 (text; ASCII)");

  NA_character_
})


setMethodS3("as.character", "AffymetrixCdfFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("File format: %s", getFileFormat(this)));
  s <- c(s, sprintf("Dimension: %s", paste(getDimension(this), collapse="x")));
  s <- c(s, sprintf("Number of cells: %d", nbrOfCells(this)));
  # Requires reading of data:
#  nbrOfPms <- sum(isPm(this));
#  s <- c(s, sprintf("Number of PM cells: %d (%.2f%%)", nbrOfPms, 100*nbrOfPms/nbrOfCells(this)));
  s <- c(s, sprintf("Number of units: %d", nbrOfUnits(this)));
  s <- c(s, sprintf("Cells per unit: %.2f", nbrOfCells(this)/nbrOfUnits(this)));
  # Requires that unit names are read:
#  s <- c(s, sprintf("Number of AFFX- units: %d", length(indexOf(this, "^AFFX-"))));
  s <- c(s, sprintf("Number of QC units: %d", nbrOfQcUnits(this)));
  s;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixCdfFile object from a CDF file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{filename}{The filename of to the file.}
#  \item{path}{The path to the file.}
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns an instance of @see "AffymetrixCdfFile" or its subclasses.
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "AffymetrixCdfFile", function(static, filename, path=NULL, ...) {
  # Arguments 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path,
                                                              mustExist=TRUE);

  # Assert that it is a CDF file
  header <- .readCdfHeader(pathname);

  NextMethod("fromFile", filename=pathname);
}, static=TRUE, protected=TRUE)



setMethodS3("getDefaultExtension", "AffymetrixCdfFile", function(static, ...) {
  "cdf";
}, static=TRUE, protected=TRUE);


###########################################################################/**
# @RdocMethod findByChipType
#
# @title "Locates a CDF file from its chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{tags}{An optional @character @vector of tags.}
#  \item{pattern}{An optional @character string.}
#  \item{...}{Not used.}
#  \item{.useAffxparser}{If @TRUE, @see "affxparser::findCdf" is used if
#    the CDF could not be located.}
# }
#
# \value{
#  Returns a pathname as a @character string to the first CDF file found.
#  If non CDF with requested chip type was found, @NULL is returned.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("findByChipType", "AffymetrixCdfFile", function(static, chipType, tags=NULL, pattern=NULL, ..., .useAffxparser=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  if (is.null(chipType)) {
    # Nothing to do, e.g. may be called via findCdf()
    return(NULL);
  }
  chipType <- Arguments$getCharacter(chipType);


  args <- list(pattern=pattern);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filename extension pattern to be searched for
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ext <- getDefaultExtension(static);
  extPattern <- sprintf("[.](%s|%s)", tolower(ext), toupper(ext));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Handle deprecated <chipType>-monocell CDFs specially
  # (This part of the code will not be updated anymore /HB 2007-12-08)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pattern <- "-monocell$";
  if (regexpr(pattern, chipType) != -1) {
    newChipType <- gsub("-monocell$", ",monocell", chipType);
    parentChipType <- gsub("-monocell$", "", chipType);  # Remove tags

    # First, see if there is a new monocell, then use that
    pattern <- sprintf("^%s%s$", newChipType, extPattern);
    pathname <- findAnnotationDataByChipType(parentChipType, pattern);

    # Second, see if the old-named monocell is there
    if (is.null(pathname)) {
      pattern <- sprintf("^%s%s$", chipType, extPattern);
      pathname <- findAnnotationDataByChipType(parentChipType, pattern=pattern);
      if (!is.null(pathname)) {
        msg <- paste("Deprecated filename of monocell CDF detected. Rename CDF file by replacing dash ('-') after 'monocell' with a comma (','): ", pathname, sep="");
        throw(msg);
      }
    }

    throw("Detected obsolete filename pattern. Monocell CDF should no longer be named '.*-monocell.CDF' but rather '.*,monocell.CDF': ", pattern);
  }


  # Create the fullname
  fullname <- paste(c(chipType, tags), collapse=",");

  # Extract the name and the tags
  parts <- unlist(strsplit(fullname, split=",", fixed=TRUE));
  chipType <- parts[1];
  tags <- parts[-1];

  # Fullname pattern
  args <- list(
    chipType=chipType,
    pattern=sprintf("^%s%s$", fullname, extPattern),
    ...
  );
  pathname <- do.call(findAnnotationDataByChipType, args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    args <- list(
      chipType=chipType,
      pattern=sprintf("^%s%s[.]lnk$", fullname, extPattern),
      ...
    );
    pathname <- do.call(findAnnotationDataByChipType, args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- Arguments$getReadablePathname(pathname, mustExist=FALSE);
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the header of the CDF file"
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
#  Returns a @list structure as returned by @see "affxparser::readCdfHeader".
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getHeader", "AffymetrixCdfFile", function(this, ...) {
  header <- this$.header
  if (is.null(header)) {
    if (!isFile(this)) return(NULL)
    header <- .readCdfHeader(getPathname(this))
    this$.header <- header
  }
  header
}, private=TRUE)


setMethodS3("getPlatform", "AffymetrixCdfFile", function(this, ...) {
  "Affymetrix";
})


setMethodS3("getChipType", "AffymetrixCdfFile", function(this, fullname=TRUE, ...) {
  if (!isFile(this)) return(NA_character_)
  chipType <- getHeader(this)$chiptype;

  # Get the main chip type?
  if (!fullname) {
    # Handle '-monocell' specially
    pattern <- "^(.*)-(monocell)$";
    if (regexpr(pattern, chipType) != -1) {
      chipType <- gsub(pattern, "\\1", chipType);
      tags <- "monocell";
    } else {
      name <- gsub("[,].*$", "", chipType);

      # Keep anything after the data-set name (and the separator).
      tags <- substring(chipType, nchar(name)+2);
      tags <- unlist(strsplit(tags, split=",", fixed=TRUE));
      if (length(tags) == 0)
        tags <- NULL;

      chipType <- name;
    }
    attr(chipType, "tags") <- tags;
  }

  chipType;
})

setMethodS3("getDimension", "AffymetrixCdfFile", function(this, ...) {
  if (!isFile(this)) return(c(NA_integer_, NA_integer_))
  header <- getHeader(this);
  c(nbrOfRows=header$rows, nbrOfColumns=header$cols);
})

setMethodS3("nbrOfRows", "AffymetrixCdfFile", function(this, ...) {
  as.integer(getDimension(this, ...)[1]);
})

setMethodS3("nbrOfColumns", "AffymetrixCdfFile", function(this, ...) {
  as.integer(getDimension(this, ...)[2]);
})


###########################################################################/**
# @RdocMethod getUnitNames
#
# @title "Gets the names of each unit"
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
# }
#
# \value{
#  Returns a @vector of @character names.
# }
#
# \details{
#   Once read from file, this information is cached in memory for efficiency.
#   The cache can be cleared by calling \code{gc(cdf)}.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getUnitNames", "AffymetrixCdfFile", function(this, units=NULL, ...) {
  names <- this$.unitNames;

  if (is.null(names)) {
    names <- .readCdfUnitNames(getPathname(this), ...);
    this$.unitNames <- names;
  }

  if (!is.null(units))
    names <- names[units];

  names;
})


setMethodS3("hasUnitTypes", "AffymetrixCdfFile", function(this, types, ..., verbose=FALSE) {
  # Argument 'types':
  types <- Arguments$getIntegers(types, range=c(0,99));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Already have unit types cached?
  unitTypes <- this$.unitTypes;
  if (!is.null(unitTypes)) {
    hasUnitTypes <- any(unitTypes %in% types);
    return(hasUnitTypes);
  }

  # ...otherwise, scan for unit types
  allUnits <- seq_len(nbrOfUnits(this));
  chunkSize <- 5000;
  while (length(allUnits) > 0) {
    idxs <- 1:min(chunkSize, length(allUnits));
    units <- allUnits[idxs];
    unitTypes <- getUnitTypes(this, units=units, .cache=FALSE);
    hasUnitTypes <- any(unitTypes %in% types);
    if (hasUnitTypes)
      return(TRUE);
    allUnits <- allUnits[-idxs];
  }

  FALSE;
})



###########################################################################/**
# @RdocMethod getUnitTypes
#
# @title "Gets the types of a set of units"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{map}{A @character string specifying the mapping used.}
#   \item{...}{Not used.}
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
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getUnitTypes", "AffymetrixCdfFile", function(this, units=NULL, ..., force=FALSE, .cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asUnitTypeIndex <- function(unitType) {
    # From the Fusion SDK documentation:
    # ASCII:
    #   0 - Unknown, 1 - CustomSeq, 2 - Genotyping, 3 - Expression,
    #   7 - Tag/GenFlex, 8 - Copy Number
    # XDA/binary:
    #   1 - Expression, 2 - Genotyping, 3 - CustomSeq, 4 - Tag,
    #   5 - Copy Number

    map <- c("unknown"=0, "expression"=1, "genotyping"=2, "resequencing"=3, "tag"=4, "copynumber"=5, "genotypingcontrol"=6, "expressioncontrol"=7);
    storage.mode(map) <- "integer";

    res <- match(unitType, names(map)) - as.integer(1);
    attr(res, "typeMap") <- map;

    res;
  } # asUnitTypeIndex()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  types <- this$.unitTypes;


  if (.cache) {
    if (force || is.null(types)) {
      # Check in file cache
      chipType <- getChipType(this);
      key <- list(method="getUnitTypes", class=class(this)[1], version="2008-09-03", chipType=chipType);
      if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
        key <- getCacheKey(this, method="getUnitTypes", chipType=chipType);
      }
      dirs <- c("aroma.affymetrix", chipType);
      if (force) {
        types <- NULL;
      } else {
        types <- loadCache(key=key, dirs=dirs);
      }

      if (is.null(types)) {
        verbose && enter(verbose, "Reading types for *all* units");

## ISSUE: readCdfUnits() does not translate the unit types, which means
## that the unit type integer different for ASCII and binary CDFs.
## types <- readCdfUnits(getPathname(this), readType=TRUE, readDirection=FALSE, readIndices=FALSE, readXY=FALSE, readBases=FALSE, readExpos=FALSE);
## WORKAROUND: Use readCdf() which return unit type strings.
## Requires: affxparser v1.13.5 or newer.

        types <- .readCdf(getPathname(this), readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE, readUnitType=TRUE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=FALSE, readIsPm=FALSE);
        types <- unlist(types, use.names=FALSE);

        # Sanity check
        if (length(types) != nbrOfUnits(this)) {
          throw("Internal error: Number of read unit types does not match the number of units in the CDF: ", length(types), " != ", nbrOfUnits(this));
        }

        # Translate
        types <- asUnitTypeIndex(types);

        saveCache(types, key=key, dirs=dirs);

        this$.unitTypes <- types;

        verbose && exit(verbose);
      }
    }

    if (!is.null(units)) {
      map <- attr(types, "typeMap");
      types <- types[units];
      attr(types, "typeMap") <- map;
    }
  } else {
## ISSUE: types <- readCdfUnits(getPathname(this), units=units, readType=TRUE, readDirection=FALSE, readIndices=FALSE, readXY=FALSE, readBases=FALSE, readExpos=FALSE);
    types <- .readCdf(getPathname(this), readXY=FALSE, readBases=FALSE, readIndexpos=FALSE, readAtoms=FALSE, readUnitType=TRUE, readUnitDirection=FALSE, readUnitNumber=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE, readGroupDirection=FALSE, readIndices=FALSE, readIsPm=FALSE);
    types <- unlist(types, use.names=FALSE);
    types <- asUnitTypeIndex(types);
  }

  typeMap <- attr(types, "typeMap");
  types <- as.integer(types);
  attr(types, "types") <- typeMap;

  types;
})



setMethodS3("getGroupDirections", "AffymetrixCdfFile", function(this, units=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  groupDirections <- this$.groupDirections;

  if (force || is.null(groupDirections)) {
    # Check in file cache
    chipType <- getChipType(this);
    key <- list(method="getGroupDirections", class=class(this)[1], chipType=chipType);
    if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
      key <- getCacheKey(this, method="getGroupDirections", chipType=chipType);
    }
    dirs <- c("aroma.affymetrix", chipType);
    if (force) {
      groupDirections <- NULL;
    } else {
      groupDirections <- loadCache(key=key, dirs=dirs);
    }

    if (is.null(groupDirections)) {
      verbose && enter(verbose, "Reading directions for *all* unit groups");
      # Have to read some group field in order to get group directions
      groupDirections <- .readCdfUnits(getPathname(this), readExpos=TRUE,
        readBases=FALSE, readXY=FALSE, readType=FALSE, readDirection=TRUE);

      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      # Remove all 'expos'
      verbose && enter(verbose, "Removing all 'expos'");
      groupDirections <- lapply(groupDirections, FUN=function(unit) {
        groups <- .subset2(unit, 2);
        groups <- lapply(groups, FUN=.subset, 2);
        list(groups=groups);
      });

      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      verbose && enter(verbose, "Restructuring");
      groupDirections <- restruct(this, groupDirections,
                                                 verbose=less(verbose, 5));

      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      verbose && enter(verbose, "Unlisting each unit");
      groupDirections <- lapply(groupDirections, FUN=unlist,
                                                          use.names=FALSE);

      gc <- gc();
      verbose && print(verbose, gc);
      verbose && exit(verbose);

      saveCache(groupDirections, key=key, dirs=dirs);
      verbose && exit(verbose);
    }

    this$.groupDirections <- groupDirections;
  }

  if (!is.null(units))
    groupDirections <- groupDirections[units];

  groupDirections;
}, private=TRUE)



###########################################################################/**
# @RdocMethod getCellIndices
#
# @title "Gets the cell indices unit by unit"
#
# \description{
#  @get "title" of all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Additional arguments passed to
#      @see "affxparser::readCdfCellIndices".}
#   \item{useNames}{If @TRUE, element names are returned, otherwise not.}
#   \item{unlist}{If @TRUE, the unlisted result is returned. Using this
#      argument is more memory efficient that calling @see "base::unlist"
#      afterwards.}
#   \item{force}{If @TRUE, cached values are ignored.}
#   \item{cache}{If @TRUE, results are cached, if not too large.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @list structure returned by
#  @see "affxparser::readCdfCellIndices".
# }
#
# \details{
#  Note, that it is much more memory efficient to do
#  \code{getCellIndices(cdf, useNames=FALSE, unlist=TRUE)}
#  compare with \code{unlist(getCellIndices(cdf), use.names=FALSE)}.
# }
#
# \seealso{
#   See @seemethod "setRestructor" to set a default re-constructor for
#   the returned CDF structure.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCellIndices", "AffymetrixCdfFile", function(this, units=NULL, ..., useNames=TRUE, unlist=FALSE, force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getCellIndicesChunk <- function(pathname, ..., verbose=FALSE) {
    verbose && enter(verbose, "Querying CDF file");
    cdfChunk <- .readCdfCellIndices(pathname, ...);
    verbose && exit(verbose);

    # Garbage collect
    gc <- gc();

    verbose && enter(verbose, "Restructuring");
    # Always call restruct() after a readCdfNnn()!
    cdfChunk <- restruct(this, cdfChunk, verbose=less(verbose, 5));
    verbose && exit(verbose);

    gc <- gc();
    verbose && print(verbose, gc);

    cdfChunk;
  } # getCellIndicesChunk()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, range=c(1, nbrOfCells(this)));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="getCellIndices", class=class(this)[1],
             chipType=getChipType(this), units=units, ...,
             useNames=useNames, unlist=unlist);
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getCellIndices", chipType=getChipType(this), units=units, ..., useNames=useNames, unlist=unlist);
  }

  # This is a trick to store either to memory or file cache
  key <- getChecksum(key);
  if (!force) {
    # (a) Check memory cache
    res <- this$.cellIndices[[key]];

##    # (b) Check file cache
##    if (is.null(res)) {
##      res <- loadCache(key=list(key));
##    }

    if (!is.null(res)) {
      verbose && cat(verbose, "getCellIndices.AffymetrixCdfFile(): Returning cached data");
      return(res);
    }
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read from CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading cell indices from CDF file");
  verbose && cat(verbose, "Pathname: ", getPathname(this));
  verbose && cat(verbose, "Units: ");
  verbose && str(verbose, units);
  verbose2 <- -as.integer(verbose)-1;

  units0 <- units;
  if (is.null(units)) {
    units <- seq_len(nbrOfUnits(this));
  }
  nbrOfUnits <- length(units);


  if (unlist) {
    cdf <- lapplyInChunks(units, function(unitsChunk) {
      cdfChunk <- getCellIndicesChunk(getPathname(this), units=unitsChunk, ..., verbose=verbose);
      res <- vector("list", length(unitsChunk));
      res[[1]] <- unlist(cdfChunk, use.names=useNames);
      res;
    }, chunkSize=100e3, useNames=useNames, verbose=verbose);
    cdf <- unlist(cdf, use.names=useNames);
  } else {
    cdf <- lapplyInChunks(units, function(unitsChunk) {
      getCellIndicesChunk(getPathname(this), units=unitsChunk, ..., verbose=verbose);
    }, chunkSize=100e3, useNames=useNames, verbose=verbose);
  }


  verbose && exit(verbose);

  units <- units0;
  # Not needed anymore
  units0 <- NULL;

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Store read units in cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    verbose && cat(verbose, "readUnits.AffymetrixCdfFile(): Updating cache");
    # Cache small objects in memory
    if (object.size(cdf) < 10e6) {
      this$.cellIndices <- list();
      this$.cellIndices[[key]] <- cdf;
##    } else {
##      saveCache(cdf, key=list(key));
    }
  }

  cdf;
}) # getCellIndices()



setMethodS3("restruct", "AffymetrixCdfFile", function(this, cdf, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Rearrange CDF structure?
  fcn <- this$.restructor;
  if (!is.null(fcn)) {
    verbose && enter(verbose, "Restructuring CDF list structure");
    verbose && cat(verbose, "Restructuring function:");
    verbose && str(verbose, fcn);
    verbose && cat(verbose, "First element in list to be restructured:");
    verbose && str(verbose, cdf[1]);

    cdf <- fcn(cdf);

    verbose && exit(verbose);
  }

  cdf;
}, private=TRUE)





###########################################################################/**
# @RdocMethod setRestructor
#
# @title "Specifies a function through which"
#
# \description{
#  @get "title" of all or a subset of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units of interest. If @NULL, all units are considered.}
#   \item{...}{Additional arguments passed to
#      @see "affxparser::readCdfCellIndices".}
# }
#
# \value{
#  Returns the @list structure returned by
#  @see "affxparser::readCdfCellIndices".
# }
#
# \section{Requirements}{
#   The reconstructor function \emph{have to}:
#  \enumerate{
#   \item Accept a @list structure as its first argument.
#   \item Require no other arguments.
#   \item Return a @list structure of identical length as the input one.
#         In other words, it must not change the number of units nor
#         the order of units.
#  }
#
#  The reconstructor function \emph{may}:
#  \enumerate{
#   \item Rearrange the groups or change the number of groups, for
#         instance by utilizing @see "affxparser::applyCdfGroups".
#   \item Exclude some cells, for instance by excluding MMs.
#  }
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("setRestructor", "AffymetrixCdfFile", function(this, fcn=NULL, ...) {
  if (is.null(fcn)) {
  } else if (is.function(fcn)) {
  } else {
    throw("Argument 'fcn' must be NULL or a function: ", mode(fcn));
  }
  if (!identical(this$.restructor, fcn)) {
    this$.restructor <- fcn;
    clearCache(this);
  }
  invisible(this);
}, private=TRUE)

setMethodS3("getRestructor", "AffymetrixCdfFile", function(this, ...) {
  this$.restructor;
}, private=TRUE)



###########################################################################/**
# @RdocMethod readUnits
#
# @title "Reads CDF data unit by unit"
#
# \description{
#  @get "title" for all or a subset of units (probeset).
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Additional arguments passed to @see "affxparser::readCdfUnits".}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @list structure that @see "affxparser::readCdfUnits" returns
#  (possibly restructured).
# }
#
# \section{Caching}{
#   CDF data is neither cached in memory nor on file by this method.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
# NOTE: getUnits() does not work because an S4 class stole it!!!
setMethodS3("readUnits", "AffymetrixCdfFile", function(this, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  cdf <- .readCdfUnits(filename=getPathname(this), units=units, ...);

  # Always call restruct() after a readCdfNnn()!
  restruct(this, cdf, verbose=less(verbose, 5));
})




###########################################################################/**
# @RdocMethod isPm
#
# @title "Checks which cells (probes) are PMs and not"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{The units to be read. If @NULL, all units are read.}
#   \item{...}{Additional arguments passed to @see "affxparser::readCdfUnits".}
#   \item{force}{If @TRUE, cached results are ignored.}
#   \item{cache}{If @TRUE, results are cached.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @logical @vector of length K, where K equals the total number
#  of cells in the requested units.  Note that the cells are ordered as
#  they occur in the units, that is, \emph{not} in incremental order.
# }
#
# \section{Caching}{
#   This method caches a @logical @vector of length N, when N equals the
#   number of cells on the array. The size of this vector is approximately
#   4*N bytes.  The vector indicates if a cell is a PM or not.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("isPm", "AffymetrixCdfFile", function(this, units=NULL, force=FALSE, cache=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  isPm <- this$.isPm;
  if (force || is.null(isPm)) {
    if (cache) {
      # If caching, read all units
      cdf <- .readCdfIsPm(getPathname(this));
      # Always call restruct() after a readCdfNnn()!
      cdf <- restruct(this, cdf, verbose=less(verbose, 5));
      isPm <- this$.isPm <- cdf;
    } else {
      # ...otherwise, read only a subset of units
      cdf <- .readCdfIsPm(getPathname(this), units=units);
      # Always call restruct() after a readCdfNnn()!
      cdf <- restruct(this, cdf, verbose=less(verbose, 5));
      isPm <- cdf;
    }
  }

  if (cache && !is.null(units)) {
    isPm <- isPm[units];
  }

  # Return a vector
  isPm <- unlist(isPm, use.names=FALSE);

  isPm;
})


setMethodS3("identifyCells", "AffymetrixCdfFile", function(this, indices=NULL, from=1, to=nbrOfCells(this), types=c("all", "pmmm", "pm", "mm", "qc"), ..., sort=TRUE, .force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'types':
  if (is.null(types))
    types <- "all";

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Arguments 'from' and 'to':
  nbrOfCells <- nbrOfCells(this);
  from <- Arguments$getInteger(from, range=c(1, nbrOfCells));
  to <- Arguments$getInteger(to, range=c(1, nbrOfCells));

  # Argument 'indices':
  if (is.numeric(indices)) {
    getFraction <- (length(indices) == 1 && indices > 0 && indices < 1);
    if (getFraction) {
      by <- 1/indices;
    } else {
      indices <- Arguments$getIntegers(indices, range=c(1, nbrOfCells));
    }
  } else {
    getFraction <- FALSE;
  }

  if ("all" %in% types) {
    other <- 1:nbrOfCells;
  } else {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Check for cached results (already here)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Create a cache key (already here)
    verbose && enter(verbose, "Checking cache");
    chipType <- getChipType(this);
    key <- list(method="identifyCells", class=class(this)[1], chipType=chipType, indices=indices, from=from, to=to, types=types, sort=sort);
    if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
      key <- getCacheKey(this, method="identifyCells", chipType=chipType, indices=indices, from=from, to=to, types=types, sort=sort);
    }
    comment <- sprintf("%s: %s", key$method, key$chipType);
    dirs <- c("aroma.affymetrix", chipType);
    if (!.force) {
      cache <- loadCache(key=key, dirs=dirs);
      if (!is.null(cache)) {
        verbose && exit(verbose, suffix="...hit");
        return(cache);
      }
    }
    verbose && exit(verbose);
  }

  # Argument 'from':
  if (is.null(indices)) {
    indices <- seq(from=from, to=to, ...);
    indices <- as.integer(indices+0.5);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Intersect 'indices' and 'types'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!"all" %in% types) {
    verbose && enter(verbose, "Identifies cells of certain kind");
    verbose && cat(verbose, "Indices:");
    verbose && str(verbose, indices);

    indices <- getCellIndices(this, useNames=FALSE, unlist=TRUE,
                                            verbose=less(verbose));

    other <- c();
    for (type in types) {
      if (type == "pm") {
        verbose && cat(verbose, "Using PM only");
        other <- c(other, indices[isPm(this)]);
      } else if (type == "mm") {
        verbose && cat(verbose, "Using MM only");
        other <- c(other, indices[!isPm(this)]);
      } else if (type == "pmmm") {
        verbose && cat(verbose, "Using PM & MM");
        other <- c(other, indices);
      } else if (type == "qc") {
        verbose && cat(verbose, "Using QC cells only");
        # Get cell indices for all non-regular units, i.e. QCs
        other <- c(other, setdiff(1:nbrOfCells, indices));
      }
    }

    other <- unique(other);
    verbose && exit(verbose);
  } # if (!"all" ...)

  if (is.null(indices)) {
    indices <- other;
    # Not needed anymore
    other <- NULL;
  } else {
    if (getFraction) {
      # Get the fraction from the already filtered cell indices
      indices <- other[seq(from=1, to=length(other), by=by)];
    } else if (!"all" %in% types) {
      indices <- intersect(indices, other);
    }
  }

  if (sort)
    indices <- sort(indices);

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save result to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!"all" %in% types) {
    saveCache(indices, key=key, comment=comment, dirs=dirs);
  }

  indices;
}, private=TRUE);


setMethodS3("getFirstCellIndices", "AffymetrixCdfFile", function(this, units=NULL, stratifyBy=NULL, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Trying to load cached results");
  chipType <- getChipType(this);
  key <- list(method="getFirstCellIndices", class=class(this)[1], chipType=chipType, stratifyBy=stratifyBy, restructor=body(this$.restructor));
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getFirstCellIndices", chipType=chipType, stratifyBy=stratifyBy, restructor=body(this$.restructor));
  }
  dirs <- c("aroma.affymetrix", chipType);
  res <- if (force) {
    NULL;
  } else {
    loadCache(key=key, dirs=dirs);
  }
  verbose && exit(verbose);

  if (is.null(res)) {
    verbose && enter(verbose, "Reading all cell indices (slow)");
    res <- getCellIndices(this, units=NULL, ..., stratifyBy=stratifyBy, verbose=verbose);
    verbose && exit(verbose);

    verbose && enter(verbose, "Extracting the first cell in each unit group");
    # For each unit and each group, get the index of the first cell.
    res <- .applyCdfGroups(res, function(groups) {
      # For each group, pull out the first cell.
      lapply(groups, FUN=function(group) {
        # group$indices[1] == group[[1]][1] == ...
        list(indices=.subset(.subset2(group, 1), 1));
      })
    });
    verbose && exit(verbose);

    # Save to cache file
    verbose && enter(verbose, "Saving results to cache");
    saveCache(res, key=key, dirs=dirs);
    verbose && exit(verbose);
  }

  # Subset?
  if (!is.null(units))
    res <- res[units];

  res;
}, private=TRUE)


###########################################################################/**
# @RdocMethod compare
#
# @title "Checks if two AffymetrixCdfFile objects are equal"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{other}{The other @see "AffymetrixCdfFile" object to be compared to.}
#   \item{...}{Additional arguments passed to @see "affxparser::compareCdfs".}
# }
#
# \value{
#  Returns @TRUE if the two objects are equal, otherwise @FALSE.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("compare", "AffymetrixCdfFile", function(this, other, ...) {
  if (!inherits(other, "AffymetrixCdfFile"))
    return(FALSE);

  # Check if it is the same object
  if (equals(this, other))
    return(TRUE);

  res <- .compareCdfs(getPathname(this), getPathname(other), ...);

  res;
}, private=TRUE)



###########################################################################/**
# @RdocMethod convert
#
# @title "Converts a CDF into the same CDF but with another format"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{chipType}{The chip type of the new CDF.}
#   \item{suffix}{A suffix added to the chip type of the new CDF.}
#   \item{sep}{A string separating the chip type and the suffix string.}
#   \item{path}{The path where to store the new CDF file.}
#   \item{...}{Additional arguments passed to @see "affxparser::convertCdf".}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the new CDF as an @see "AffymetrixCdfFile" object.
# }
#
# \seealso{
#   To compare two CDFs, use \code{equals()}.
#   Internally @see "affxparser::convertCdf" is used.
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("convert", "AffymetrixCdfFile", function(this, chipType=getChipType(this), suffix=NULL, sep="-", path="cdf", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Create the pathname of the destination CDF
  name <- paste(c(chipType, suffix), collapse=sep);
  dest <- sprintf("%s.CDF", name);
  dest <- Arguments$getWritablePathname(dest, path=path);

  # Convert CDF
  src <- getPathname(this);
  verbose2 <- -getThreshold(verbose);
  res <- .convertCdf(src, dest, ..., verbose=verbose2);

  # Return an AffymetrixCdfFile object for the new CDF
  cdf <- newInstance(this, dest)

  ## Create checksum
  cdfZ <- getChecksumFile(cdf)

  cdf
})




###########################################################################/**
# @RdocMethod getGenomeInformation
#
# @title "Gets genome information for this chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{types}{A @character @vector specifying what type of genome
#     information sets to search for.}
#   \item{...}{Not used.}
#   \item{force}{If @FALSE, cached information is retrieved, otherwise not.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "GenomeInformation" object.
# }
#
# \examples{\dontrun{
#   @include "../incl/getGenomeInformation.Rex"
# }}
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getGenomeInformation", "AffymetrixCdfFile", function(this, types=c("UGP", "dChip"), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'types':
  types <- Arguments$getCharacters(types);
  types <- tolower(types);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating a GenomeInformation file");

  chipType <- getChipType(this, fullname=FALSE);
  tags <- getTags(this);
  tags <- setdiff(tags, "monocell");
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=", "));
  nbrOfUnits <- nbrOfUnits(this);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);

  gi <- this$.gi;
  if (force || is.null(gi)) {
    gi <- NULL;
    for (type in types) {
      verbose && enter(verbose, "Searching for type '", type, "'");

      tryCatch({
        if (type == "ugp") {
          gi <- UgpGenomeInformation$byChipType(chipType, tags=tags,
                     nbrOfUnits=nbrOfUnits, verbose=less(verbose, 5));
          break;
        } else if (type == "dchip") {
          gi <- DChipGenomeInformation$byChipType(chipType,
                     nbrOfUnits=nbrOfUnits, verbose=less(verbose, 5));

          break;
        }
      }, error = function(ex) {})

      verbose && exit(verbose);
    } # for (type ...)

    if (is.null(gi)) {
      throw("Failed to retrieve genome information for this chip type: ",
                                                               chipType);
    }

    setGenomeInformation(this, gi);
  }

  verbose && exit(verbose);

  gi;
})


setMethodS3("setGenomeInformation", "AffymetrixCdfFile", function(this, gi=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'gi':
  if (!is.null(gi)) {
    gi <- Arguments$getInstanceOf(gi, "GenomeInformation");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isCompatible <- isCompatibleWithCdf(gi, this);
  if (!isCompatible) {
    throw("Cannot set genome information. The object 'gi' ('", getFullName(gi), "') of ", class(gi)[1], " is not compatible with the CDF ('", getFullName(this), "'). The reason was: ", attr(isCompatible, "reason"));
  }


  this$.gi <- gi;

  invisible(this);
}, protected=TRUE)


setMethodS3("getSnpInformation", "AffymetrixCdfFile", function(this, types=c("UFL", "dChip"), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'types':
  types <- Arguments$getCharacters(types);
  types <- tolower(types);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating a SnpInformation file");

  chipType <- getChipType(this, fullname=FALSE);
  tags <- getTags(this);
  tags <- setdiff(tags, "monocell");
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Tags: ", paste(tags, collapse=", "));
  nbrOfUnits <- nbrOfUnits(this);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);

  si <- this$.si;
  if (force || is.null(si)) {
    gi <- NULL;
    for (type in types) {
      verbose && enter(verbose, "Searching for type '", type, "'");

      tryCatch({
        if (type == "ufl") {
          si <- UflSnpInformation$byChipType(chipType, tags=tags,
                     nbrOfUnits=nbrOfUnits, verbose=less(verbose, 5));
          break;
        } else if (type == "dchip") {
          si <- DChipSnpInformation$byChipType(chipType,
                     nbrOfUnits=nbrOfUnits, verbose=less(verbose, 5));
          break;
        }
      }, error = function(ex) {})

      verbose && exit(verbose);
    } # for (type ...)

    if (is.null(si)) {
      throw("Failed to retrieve SNP information for this chip type: ",
                                                            chipType);
    }

    setSnpInformation(this, si);
  }

  verbose && exit(verbose);

  si;
}, private=TRUE)


setMethodS3("setSnpInformation", "AffymetrixCdfFile", function(this, si=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'si':
  if (!is.null(si)) {
    si <- Arguments$getInstanceOf(si, "SnpInformation");
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Sanity checks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  isCompatible <- isCompatibleWithCdf(si, this);
  if (!isCompatible) {
    throw("Cannot set genome information. The object 'si' ('", getFullName(si), "') of ", class(si)[1], " is not compatible with the CDF ('", getFullName(this), "'). The reason was: ", attr(isCompatible, "reason"));
  }

  this$.si <- si;

  invisible(this);
}, protected=TRUE)


###########################################################################/**
# @RdocMethod convertUnits
#
# @title "Gets and validates unit indices"
#
# \description{
#  @get "title" either by unit names or by a unit indices (validation).
# }
#
# @synopsis
#
# \arguments{
#   \item{units}{Either a @character @vector with unit names, or an @integer
#     @vector with unit indices to be validated.
#     If @NULL, all unit indices are returned.}
#   \item{keepNULL}{If @TRUE, @NULL returns @NULL.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector with unit indices.
#  If some units are non existing, an error is thrown.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("convertUnits", "AffymetrixCdfFile", function(this, units=NULL, keepNULL=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'units':
  if (is.null(units)) {
    # Return all units
    if (keepNULL)
      return(NULL);
    units <- 1:nbrOfUnits(this);
  } else if (is.character(units)) {
    # Identify units by their names
    unitNames <- units;
    units <- indexOf(this, names=unitNames);
    missing <- unitNames[is.na(units)];
    n <- length(missing);
    if (n > 0) {
      throw(sprintf("Argument 'units' contains unknown unit names: %s [%d]",
                                                      hpaste(missing), n));
    }
  } else {
    # Validate unit indices
    units <- Arguments$getIndices(units, max=nbrOfUnits(this));
  }

  units;
}, private=TRUE)


setMethodS3("validate", "AffymetrixCdfFile", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  assertUnits <- function(expr, fmtstr="%d unit(s) (i.e. %s) are invalid: %s") {
    units <- which(expr);
    nunits <- length(units);
    if (nunits > 0L) {
      fmtstr <- paste("Detected invalid/corrupt CDF: ", fmtstr, sep="");
      msg <- sprintf(fmtstr, nunits, hpaste(units), getPathname(this));
      throw(msg);
    }
  } # assertUnits()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for units with zero unit groups
  #
  # Examples:
  # o HTHGU133A_Hs_ENTREZG.cdf (v 12.0.0):
  #    Error: Detected 1 unit(s) (i.e. 11973) with zero unit groups: ...
  #   because it's CDF header claims to have 11,973 units, whereas there
  #   are only 11,972.  See also thread '[customcdf] ENTREZG, AUGUSTUST
  #   for pig species is updated' on May 8, 2012 [http://goo.gl/Xg1pp]
  #
  # Examples:
  # o HTHGU133A_Hs_ENTREZG.cdf (v 12.0.0) [as above]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ns <- nbrOfGroupsPerUnit(this);
  assertUnits(ns == 0L, "%d unit(s) (i.e. %s) with zero unit groups: %s");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for empty unit names
  #
  # Examples:
  # o HTHGU133A_Hs_ENTREZG.cdf (v 12.0.0):
  #    Error: Detected 1 unit(s) (i.e. 11973) with empty unit names: ...
  #   because it's CDF header claims to have 11,973 units, whereas there
  #   are only 11,972.  See also thread '[customcdf] ENTREZG, AUGUSTUST
  #   for pig species is updated' on May 8, 2012 [http://goo.gl/Xg1pp]
  #
  # Examples:
  # o HTHGU133A_Hs_ENTREZG.cdf (v 12.0.0) [as above]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  unitNames <- getUnitNames(this);
  assertUnits(((ns == 0L) & (!nzchar(unitNames))), "%d unit(s) (i.e. %s) with zero unit groups and empty unit names: %s");


  invisible(this);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2012-12-16
# o Now validate() for AffymetrixCdfFile accepts empty unit names as
#   long as the unit is not empty.
# 2012-10-18
# o Added validate() for AffymetrixCdfFile, which validate a CDF for
#   the most "common" errors, to help troubleshooting.  Note that the
#   validation is not complete, i.e. rare/unknown errors are not caught.
# 2012-10-14
# o CLEANUP: findByChipType() for AffymetrixCdfFile no longer support
#   monocell CDF file named <chipType>-monocell.CDF, and gives an
#   informative error if that is still the case.  Since December 2012,
#   the filename should instead be <chipType>,monocell.CDF.
# 2011-11-18
# o CLEANUP: Now the filename extension pattern for findByChipType()
#   of AffymetrixCdfFile is inferred from getDefaultExtension().
# o Added static getDefaultExtension() for AffymetrixCdfFile.
# 2011-03-04
# o ROBUSTNESS: Now getCellIndices() for AffymetrixCdfFile asserts that
#   argument 'units' can be coerced to integer indices.
# 2010-05-09
# o Added more explicit garbage collection to getGroupDirections() for
#   the AffymetrixCdfFile class.
# 2010-01-03
# o BUG FIX: After loading aroma.affymetrix, findCdf() would give "Error in
#   if (regexpr(pattern, chipType) != -1) { : argument is of length zero",
#   because AffymetrixCdfFile$findByChipType(chipType=NULL) was not valid.
#   Now the latter returns NULL without complaining.
# 2009-07-08
# o Added getUnitTypesFile() for AffymetrixCdfFile.
# o Now AffymetrixCdfFile implements also the UnitTypesFile interface.
# 2009-05-09
# o Added names to the returned dimension of getDimension().
# 2009-02-10
# o Added selection/validation of number of units in
#   get(Genome|Snp)nformation().
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
