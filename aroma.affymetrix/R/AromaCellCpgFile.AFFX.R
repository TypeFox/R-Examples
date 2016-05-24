# @author "MR, HB"
setMethodS3("allocateFromCdf", "AromaCellCpgFile", function(static, cdf, path=NULL, tags="*", ...) {
  # Argument 'cdf':
  cdf <- Arguments$getInstanceOf(cdf, "AffymetrixCdfFile");

  # Argument 'tags':
  tags <- Arguments$getTags(tags, collapse=NULL);


  chipType <- getChipType(cdf);
  parts <- strsplit(chipType, split=",", fixed=TRUE);
  parts <- unlist(parts, use.names=FALSE);

  tags[tags == "*"] <- paste(parts[-1], collapse=",");
  chipType <- parts[1];

  # Drop unwanted tags. AD HOC. /HB 2010-02-22
  chipType <- gsub(",monocell", "", chipType, fixed=TRUE);

  # Output path
  if (is.null(path)) {
    path <- file.path("annotationData", "chipTypes", chipType);
  }
  path <- Arguments$getWritablePath(path);

  platform <- getPlatform(cdf);
  nbrOfCells <- nbrOfCells(cdf);
  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);
  allocate(static, filename=filename, path=path, nbrOfCells=nbrOfCells,
                                 platform=platform, chipType=chipType, ...);
}, static=TRUE)



############################################################################
# HISTORY:
# 2011-02-28
# o STANDARDIZATION: Now the default output path for all allocateFromCdf()
#   is annotationData/chipTypes/<chipType>/.  Before it was the same
#   directory as the CDF, which may for instance have been in a deeper
#   subdirectory, or recently also in a sibling root path.
# 2010-03-14 [HB]
# o BUG FIX: allocateFromCdf() of AromaCellCpgFile would drop all but the
#   first tag.
# 2009-02-22 [HB]
# o Forgot to make allocateFromCdf() of AromaCellCpgFile static.
# 2010-02-19 [MR]
# o Created from AromaCellPositionFile.AFFX.R.
############################################################################
