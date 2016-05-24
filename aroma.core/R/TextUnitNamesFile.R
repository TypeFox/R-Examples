setConstructorS3("TextUnitNamesFile", function(..., platform=NULL) {
  # Argument 'platform':
  if (!is.null(platform)) {
    platform <- Arguments$getCharacter(platform, length=c(1,1));
  }

  extend(TabularTextFile(...), c("TextUnitNamesFile",
    uses("UnitNamesFile"), uses("FileCacheKeyInterface")),
    .platform = platform,
    "cached:.unitNames" = NULL
  );
})


setMethodS3("as.character", "TextUnitNamesFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Number units: %d", nbrOfUnits(this, fast=TRUE)));
  s;
}, protected=TRUE)


setMethodS3("getPlatform", "TextUnitNamesFile", function(this, ..., force=FALSE) {
  platform <- this$.platform;
  if (force || is.null(platform)) {
    args <- getHeader(this)$commentArgs;
    platform <- args[["platform"]];
    this$.platform <- platform;
  }
  platform;
})

setMethodS3("getChipType", "TextUnitNamesFile", function(this, ...) {
  getName(this);
})

setMethodS3("getUnitNames", "TextUnitNamesFile", function(this, units=NULL, ...) {
  unitNames <- this$.unitNames;
  if (is.null(unitNames)) {
    # Read all unit names
    data <- readDataFrame(this, ...);
    unitNames <- data[,1,drop=TRUE];
    unitNames <- as.character(unitNames);
    this$.unitNames <- unitNames;
    this$.nbrOfUnits <- length(unitNames);
  }

  # Subsetting
  if (!is.null(units)) {
    units <- Arguments$getIndices(units, max=length(unitNames));
    unitNames <- unitNames[units];
  }

  unitNames;
})

setMethodS3("nbrOfUnits", "TextUnitNamesFile", function(this, ...) {
  nbrOfUnits <- this$.nbrOfUnits;
  if (is.null(nbrOfUnits)) {
    nbrOfUnits <- NextMethod("nbrOfUnits");
    this$.nbrOfUnits <- nbrOfUnits;
  }
  nbrOfUnits;
})

setMethodS3("getFilenameExtension", "TextUnitNamesFile", function(static, ...) {
  "txt";
}, static=TRUE, protected=TRUE);


setMethodS3("getExtensionPattern", "TextUnitNamesFile", function(static, ...) {
  "[.](txt)$";
}, static=TRUE, protected=TRUE)


setMethodS3("findByChipType", "TextUnitNamesFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get fullname, name, and tags
  fullname <- paste(c(chipType, tags), collapse=",");
  parts <- unlist(strsplit(fullname, split=","));
  # Strip 'monocell' parts
  parts <- parts[parts != "monocell"];
  chipType <- parts[1];
  tags <- parts[-1];
  fullname <- paste(c(chipType, tags), collapse=",");

  ext <- getFilenameExtension(static);
  ext <- paste(c(tolower(ext), toupper(ext)), collapse="|");
  ext <- sprintf("(%s)", ext);

  pattern <- sprintf("^%s.*,unitNames[.]%s$", fullname, ext);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;  # Override argument 'pattern'?
#  args$firstOnly <- FALSE;
#  str(args);
  pathname <- do.call(findAnnotationDataByChipType, args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- sprintf("^%s.*[.]%s[.]lnk$", chipType, ext);
    args$pattern <- pattern;
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



############################################################################
# HISTORY:
# 2012-10-16
# o CLEANUP: Dropped getHeaderParameters() for TextUnitNamesFile, which
#   is now done by getHeader() of TabularTextFile.
# 2009-05-12
# o Now TextUnitNamesFile caches all unit names in memory.
# 2009-05-11
# o BUG FIX: getPlatform() of TextUnitNamesFile would sometimes return a
#   list of length one, instead of an single character string.
# o Added getHeaderParameters().
# 2009-04-06
# o BUG FIX: getUnitNames(..., units=NULL) of TextUnitNamesFile would
#   make the object believe there are zero units in the file.
# 2009-03-23
# o Created.
############################################################################
