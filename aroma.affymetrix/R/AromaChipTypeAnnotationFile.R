###########################################################################/**
# @RdocClass AromaChipTypeAnnotationFile
#
# @title "The AromaChipTypeAnnotationFile class"
#
# \description{
#  @classhierarchy
#
#  An AromaChipTypeAnnotationFile object represents an annotation file for a
#  specific chip type.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AffymetrixFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AromaChipTypeAnnotationFile", function(...) {
  this <- extend(AffymetrixFile(...), "AromaChipTypeAnnotationFile");

  # Parse attributes (all subclasses must call this in the constructor).
  setAttributesByTags(this)

  this;
})


setMethodS3("as.character", "AromaChipTypeAnnotationFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Path: %s", getPath(this)));
  s <- c(s, sprintf("Filename: %s", getFilename(this)));

  # File size
  fileSize <- getFileSize(this, "units");
  if (!is.na(fileSize)) {
    fileSizeB <- sprintf("%.0f bytes", getFileSize(this, "numeric"));
    if (fileSizeB != fileSize) {
      fileSize <- sprintf("%s (%s)", fileSize, fileSizeB);
    }
  }
  s <- c(s, sprintf("File size: %s", fileSize));

  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  GenericSummary(s);
}, protected=TRUE)



###########################################################################/**
# @RdocMethod fromFile
#
# @title "Sets up an AromaChipTypeAnnotationFile"
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
#  Returns an instance of @see "AromaChipTypeAnnotationFile" or its subclasses.
#  If the file is not found or if it is of the wrong file format, an
#  error is thrown.
# }
#
# \seealso{
#   @seemethod "byChipType".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("fromFile", "AromaChipTypeAnnotationFile", function(static, filename, path=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' and 'path':
  pathname <- Arguments$getReadablePathname(filename, path=path,
                                                              mustExist=TRUE);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Try to define an instance of a subclass traversing bottom up.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  clazz <- Class$forName(class(static)[1]);
  for (className in rev(getKnownSubclasses(clazz))) {
    clazz <- Class$forName(className);
    tryCatch({
      res <- newInstance(clazz, pathname);
      return(res);
    }, error = function(ex) {})
  }

  newInstance(static, pathname);
}, static=TRUE, protected=TRUE)




###########################################################################/**
# @RdocMethod byChipType
# @aliasmethod byName
#
# @title "Defines an AromaChipTypeAnnotationFile object by chip type"
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
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "AromaChipTypeAnnotationFile" object.
# }
#
# \seealso{
#   @seemethod "fromFile".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("byChipType", "AromaChipTypeAnnotationFile", function(static, chipType, tags=NULL, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Setting up ", class(static)[1L], " by chip type");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search for a matching annotation file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findByChipType(static, chipType, tags=tags, verbose=less(verbose, 50));
  verbose && cat(verbose, "Located file:");

  verbose && print(verbose, pathname);

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

  res <- fromFile(static, filename=pathname, path=NULL, ...);
  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}, static=TRUE)


setMethodS3("byName", "AromaChipTypeAnnotationFile", function(static, ...) {
  byChipType(static, ...);
}, static=TRUE, protected=TRUE)




setMethodS3("getDefaultExtension", "AromaChipTypeAnnotationFile", function(static, ...) {
  # Guess the filename extension from the class name, which might be wrong
  className <- class(static)[1];

  ext <- gsub("File$", "", className);
  ext <- strsplit(ext, split="", fixed=TRUE)[[1]];
  n <- length(ext);
  pos <- which(ext == toupper(ext));
  pos <- pos[length(pos)];
  ext <- ext[seq(from=pos, to=n)];
  ext <- paste(ext, collapse="");
  ext <- tolower(ext);

  attr(ext, "note") <- sprintf("this may not be the correct extension as it was guessed from the class name '%s'", className);

  ext;
}, static=TRUE, protected=TRUE);


###########################################################################/**
# @RdocMethod findByChipType
#
# @title "Locates an annotation file by its chip type"
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
#  \item{...}{Additional arguments.}
# }
#
# \value{
#  Returns the pathname (as a @character string) to the first annotation
#  chip type file found.  If no one was found, @NULL is returned.
# }
#
# \seealso{
#   @seemethod "byChipType".
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("findByChipType", "AromaChipTypeAnnotationFile", abstract=TRUE, protected=TRUE);





###########################################################################/**
# @RdocMethod getHeader
#
# @title "Gets the header of the annotation file"
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
#  Returns a @list structure.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getHeader", "AromaChipTypeAnnotationFile", abstract=TRUE);


setMethodS3("getPlatform", "AromaChipTypeAnnotationFile", abstract=TRUE);


setMethodS3("getChipType", "AromaChipTypeAnnotationFile", abstract=TRUE);



############################################################################
# HISTORY:
# 2012-11-29
# o byChipType() for AromaChipTypeAnnotationFile called global 'verbose'.
# 2012-10-18
# o Now as.character() for AromaChipTypeAnnotationFile reports the file
#   size as "24.66 MB (25853850 bytes)" and not just as "24.66 MB".
# 2012-06-18
# o FIX: Fixed a type in an error message generated by byChipType()
#   for AromaChipTypeAnnotationFile.
# 2011-11-19
# o Now byChipType() for AromaChipTypeAnnotationFile gives an error
#   message with more information on which file it failed to locate,
#   e.g. by specifying filename extension it looked for.
# o Added default getDefaultExtension() for AromaChipTypeAnnotationFile,
#   which guesses the filename extension from the class name.
# 2008-05-09
# o Created from AffymetrixCdfFile.R.
############################################################################
