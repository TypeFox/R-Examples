###########################################################################/**
# @RdocClass AromaMicroarrayTabularBinaryFile
#
# @title "The AromaMicroarrayTabularBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaMicroarrayTabularBinaryFile is an abstract
#  @see "AromaTabularBinaryFile".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#
# \seealso{
#   @see "AromaTabularBinaryFile".
# }
#*/###########################################################################
setConstructorS3("AromaMicroarrayTabularBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), c("AromaMicroarrayTabularBinaryFile",
                                              uses("AromaPlatformInterface"))
  );
}, abstract=TRUE)



setMethodS3("as.character", "AromaMicroarrayTabularBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  s <- c(s, sprintf("Platform: %s", getPlatform(this)));
  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  n <- length(s);
  s <- s[c(1:(n-2), n, n-1)];  # TODO: '[' drops the class.
  s <- GenericSummary(s);
  s;
}, protected=TRUE)


setMethodS3("getFilenameExtension", "AromaMicroarrayTabularBinaryFile", static=TRUE, abstract=TRUE);


setMethodS3("getPlatform", "AromaMicroarrayTabularBinaryFile", function(this, ...) {
  footer <- readFooter(this);
  platform <- footer$platform;

  if (is.null(platform)) {
    # AD HOC: If there is no platform information in the file, then assume
    # it is an Affymetrix platform.  Newer files are allocated with
    # platform information.  Older files have only be created for the
    # Affymetrix platform. /HB 2008-05-18.
    platform <- "Affymetrix";
    warning(sprintf("%s does not have a 'platform' footer attribute. Assuming platform 'Affymetrix': %s", class(this)[1], getPathname(this)));
  }

  if (!is.null(platform)) {
    platform <- as.character(platform);
    platform <- unlist(strsplit(platform, split="[\t]"));
    platform <- trim(platform);
  }

  platform;
})


setMethodS3("getChipType", "AromaMicroarrayTabularBinaryFile", function(this, fullname=TRUE, .old=FALSE, ...) {
  footer <- readFooter(this);
  chipType <- footer$chipType;

  if (!missing(.old)) {
    .Deprecated("Argument '.old' is deprecated since January 2008 and will be made defunct in a future version of aroma.core.");
  }

  if (is.null(chipType)) {
    msg <- paste0("File format error: This ", class(this)[1L], " file does not contain information on chip type in the file footer.  This is because the file is of an older file format an is no longer supported.  Please update to a more recent version: ", getPathname(this));
    if (.old) {
     msg <- paste0(msg, " [Argument '.old' (== TRUE) is deprecated since January 2008 and now defunct]");
    }
    throw(msg);
  }

  if (!fullname) {
    chipType <- gsub(",.*", "", chipType);
  }
  chipType <- trim(chipType);
  if (nchar(chipType) == 0L) {
    throw("File format error: The chip type according to the file footer is empty.");
  }

  if (!is.null(chipType)) {
    chipType <- as.character(chipType);
    chipType <- unlist(strsplit(chipType, split="[\t]"));
    chipType <- trim(chipType);
  }

  chipType;
}) # getChipType()



setMethodS3("byChipType", "AromaMicroarrayTabularBinaryFile", function(static, chipType, tags=NULL, validate=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType, length=c(1,1));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Locating ", class(static)[1]);

  pathname <- findByChipType(static, chipType=chipType, tags=tags,
                                                     firstOnly=TRUE, ...);
  if (is.null(pathname)) {
    ext <- getDefaultExtension(static);
    note <- attr(ext, "note");
    msg <- sprintf("Failed to create %s object. Could not locate an annotation data file for chip type '%s'", class(static)[1], chipType);
    if (is.null(tags)) {
      msg <- sprintf("%s (without requiring any tags)", msg);
    } else {
      msg <- sprintf("%s with tags '%s'", msg, paste(tags, collapse=","));
    }
    msg <- sprintf("%s and with filename extension '%s'", msg, ext);
    if (!is.null(note)) {
      msg <- sprintf("%s (%s)", msg, note);
    }
    msg <- sprintf("%s.", msg);
    throw(msg);
  }

  verbose && cat(verbose, "Located file: ", pathname);

  # Create object
  res <- newInstance(static, pathname);

  verbose && exit(verbose);

  res;
}, static=TRUE)



setMethodS3("findByChipType", "AromaMicroarrayTabularBinaryFile", function(static, chipType, tags=NULL, ...) {
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

  pattern <- sprintf("^%s.*[.]%s$", fullname, ext);
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





setMethodS3("allocate", "AromaMicroarrayTabularBinaryFile", function(static, platform=footer$platform, chipType=footer$chipType, footer=list(), ...) {
  res <- NextMethod("allocate");

  # Write attributes to footer
  attrs <- list(platform=platform, chipType=chipType);
  for (key in names(attrs)) {
    footer[[key]] <- attrs[[key]];
  }
  writeFooter(res, footer);

  res;
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2014-06-28
# o CLEANUP: Deprecated argument '.old' of getChipType() for
#   AromaMicroarrayTabularBinaryFile.
# 2008-07-09
# o Extracted from AromaUnitTabularBinaryFile.R.
############################################################################
