###########################################################################/**
# @RdocClass AffymetrixPgfFile
#
# @title "The AffymetrixPgfFile class"
#
# \description{
#  @classhierarchy
#
#  An AffymetrixPgfFile object represents a generic Affymetrix Probe Group
#  File (PGF).
#  A PGF file "provides information about what probes are contained
#  within a probeset and information about the nature of the probes
#  necessary for analysis. The current PGF file format (version 1)
#  is only specified for expression style probesets." [1]
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaChipTypeAnnotationFile".}
# }
#
# \references{
#  [1] ... \cr
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AffymetrixPgfFile", function(...) {
  this <- extend(AromaChipTypeAnnotationFile(...), c("AffymetrixPgfFile",
            uses("UnitNamesFile", "UnitTypesFile", "AromaPlatformInterface")),
    "cached:.header" = NULL,
    "cached:.data" = NULL
  );

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("getExtensionPattern", "AffymetrixPgfFile", function(static, ...) {
  "[.](pgf|PGF)$";
}, static=TRUE, protected=TRUE)


setMethodS3("getUnitNamesFile", "AffymetrixPgfFile", function(this, ...) {
  this;
}, protected=TRUE)



setMethodS3("as.character", "AffymetrixPgfFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Dimension: %s", paste(getDimension(this), collapse="x")));
  s <- c(s, sprintf("Number of cells: %d", nbrOfCells(this)));
  s <- c(s, sprintf("Number of units: %d", nbrOfUnits(this)));
  s <- c(s, sprintf("Cells per unit: %.2f", nbrOfCells(this)/nbrOfUnits(this)));
  s;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod fromFile
#
# @title "Defines an AffymetrixPgfFile object from a PGF file"
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
#  Returns an instance of @see "AffymetrixPgfFile" or its subclasses.
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
setMethodS3("fromFile", "AffymetrixPgfFile", function(static, filename, path=NULL, ...) {
  # Arguments 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path,
                                                              mustExist=TRUE);

  # Assert that it is a PGF file
  header <- .readPgfHeader(pathname);

  NextMethod("fromFile", filename=pathname);
}, static=TRUE, protected=TRUE)



setMethodS3("getDefaultExtension", "AffymetrixPgfFile", function(static, ...) {
  "pgf";
}, static=TRUE, protected=TRUE);


###########################################################################/**
# @RdocMethod findByChipType
#
# @title "Locates a PGF file from its chip type"
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
# }
#
# \value{
#  Returns a pathname as a @character string to the first PGF file found.
#  If non PGF with requested chip type was found, @NULL is returned.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("findByChipType", "AffymetrixPgfFile", function(static, chipType, tags=NULL, pattern=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  if (is.null(chipType)) {
    # Nothing to do, e.g. may be called via findPgf()
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
# @title "Gets the header of the PGF file"
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
#  Returns a @list structure as returned by @see "affxparser::readPgfHeader".
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getHeader", "AffymetrixPgfFile", function(this, ...) {
  hdr <- this$.header;
  if (is.null(hdr)) {
    pathname <- getPathname(this);
    hdr <- .readPgfHeader(pathname)$header;
    for (ff in c("rows", "cols", "probesets", "datalines")) {
      hdr[[ff]] <- Arguments$getInteger(hdr[[ff]], .name=ff);
    }
    names <- names(hdr);
    names <- gsub("_", " ", names, fixed=TRUE);
    names <- toCamelCase(names);
    names(hdr) <- names;
    this$.header <- hdr;
  }
  hdr;
}, private=TRUE)


setMethodS3("getPlatform", "AffymetrixPgfFile", function(this, ...) {
  "Affymetrix";
})


setMethodS3("getChipType", "AffymetrixPgfFile", function(this, fullname=TRUE, ...) {
  chipType <- getHeader(this)$chipType;

  # Get the main chip type?
  if (!fullname) {
    name <- gsub("[,].*$", "", chipType);

    # Keep anything after the data-set name (and the separator).
    tags <- substring(chipType, nchar(name)+2);
    tags <- unlist(strsplit(tags, split=",", fixed=TRUE));
    if (length(tags) == 0)
      tags <- NULL;

    chipType <- name;
    attr(chipType, "tags") <- tags;
  }

  chipType;
})

setMethodS3("getDimension", "AffymetrixPgfFile", function(this, ...) {
  header <- getHeader(this);
  c(nbrOfRows=header$rows, nbrOfColumns=header$cols);
})

setMethodS3("nbrOfRows", "AffymetrixPgfFile", function(this, ...) {
  as.integer(getDimension(this, ...)[1]);
})

setMethodS3("nbrOfColumns", "AffymetrixPgfFile", function(this, ...) {
  as.integer(getDimension(this, ...)[2]);
})


setMethodS3("nbrOfCells", "AffymetrixPgfFile", function(this, ...) {
  as.integer(prod(getDimension(this, ...)));
})

setMethodS3("nbrOfUnits", "AffymetrixPgfFile", function(this, ...) {
  getHeader(this)$probesets;
})


setMethodS3("readRawData", "AffymetrixPgfFile", function(this, ...) {
  rawData <- this$.rawData;
  if (is.null(rawData)) {
    pathname <- getPathname(this);
    env <- .readPgfEnv(pathname, ...);
    env$header <- NULL;

    keys <- ls(envir=env);
    fields <- grep("^probeset", keys, value=TRUE);
    data <- mget(fields, envir=env);
    names(data) <- toCamelCase(gsub("^probeset", "", names(data)));
    data$type <- as.factor(data$type);
    unitData <- as.data.frame(data, stringsAsFactors=FALSE);
    for (ff in fields) rm(list=ff, envir=env);

    keys <- ls(envir=env);
    fields <- grep("^probe", keys, value=TRUE);
    data <- mget(fields, envir=env);
    names(data) <- toCamelCase(gsub("^probe", "", names(data)));
    data$type <- as.factor(data$type);
    cellData <- as.data.frame(data, stringsAsFactors=FALSE);
    for (ff in fields) rm(list=ff, envir=env);

    keys <- ls(envir=env);
    fields <- grep("^atom", keys, value=TRUE);
    data <- mget(fields, envir=env);
    names(data) <- toCamelCase(gsub("^atom", "", names(data)));
    atomData <- as.data.frame(data, stringsAsFactors=FALSE);
    for (ff in fields) rm(list=ff, envir=env);

    rawData <- list(unitData=unitData, cellData=cellData, atomData=atomData);
    this$.rawData <- rawData;
  }

  rawData;
}, private=TRUE)



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
#   The cache can be cleared by calling \code{gc(pgf)}.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getUnitNames", "AffymetrixPgfFile", function(this, units=NULL, ...) {
  data <- readRawData(this, ...);
  names <- data$unitData$name;
  if (!is.null(units)) {
    names <- names[units];
  }
  names;
})


############################################################################
# HISTORY:
# 2012-06-14
# o Created from AffymetrixCdfFile.R.
############################################################################
