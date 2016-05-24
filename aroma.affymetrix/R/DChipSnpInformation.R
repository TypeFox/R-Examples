###########################################################################/**
# @RdocClass DChipSnpInformation
#
# @title "The DChipSnpInformation class"
#
# \description{
#  @classhierarchy
#
#  This class represents dChip genome information files, which typically
#  contains information on nucleotide sequences and fragment lengths
#  of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "SnpInformation".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   The dChip SNP information files for various chip types can be
#   downloaded from \url{http://www.hsph.harvard.edu/cli/complab/dchip/}.
#   Put each file in a
#   directory named identically as the corresponding chip type under the
#   \emph{annotations/} directory, e.g.
#   \emph{annotations/Mapping50K\_Hind240/50k hind snp info AfAm
#   june 05 hg17.xls}.
#   Note that dChip changes the filename and file format slightly between
#   chip types, but currently the @seemethod "byChipType" basically searches
#   for files with names consisting of \code{"snp info"} or
#   \code{"snp_info"}.  At least for the most common chip types, there
#   is no need to rename the files in order for this class to recognize them.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("DChipSnpInformation", function(...) {
  this <- extend(SnpInformation(...), "DChipSnpInformation")
  if (isFile(this)) verify(this)
  this
})



setMethodS3("findByChipType", "DChipSnpInformation", function(static, chipType, version=NULL, ...) {
  # Argument 'version':
  if (is.null(version))
    version <- ".*";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pattern <- sprintf("^.*( |_)snp( |_)info(| |_).*%s[.](txt|xls)$", version);
  pathname <- findAnnotationDataByChipType(chipType, pattern);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # As a backup search the "old" style
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(pathname)) {
    path <- filePath("annotations", chipType);
    path <- Arguments$getReadablePath(path, mustExist=FALSE);

    if (isDirectory(path)) {
      pathnames <- list.files(path=path, pattern=pattern, full.names=TRUE);
      nfiles <- length(pathnames);
      if (nfiles > 1) {
        pathnames <- sort(pathnames);
        warning("Found more than one matching dChip genome information file, but returning only the last one: ", paste(pathnames, collapse=", "));
        pathnames <- rev(pathnames);
        pathname <- pathnames[1];
      }
    }
  }

  pathname;
}, static=TRUE, protected=TRUE)


###########################################################################/**
# @RdocMethod byChipType
#
# @title "Defines a DChipSnpInformation object by chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{rootPath}{A @character string specifying the root path, i.e.
#    the annotation directory.}
#  \item{version}{An optional @character string specifying the version
#    string, if more than one version is available.}
#  \item{pattern}{An optional filename pattern used to locate the
#    dChip genome file.  If @NULL, a default pattern is used.}
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "DChipSnpInformation" object.
#  If no file was not found, an error is thrown.
# }
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("byChipType", "DChipSnpInformation", function(static, chipType, version=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Defining ", class(static)[1], " from chip type");
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Version: ", version);

  # Search for the genome information file
  pathname <- findByChipType(static, chipType, version=version, ..., verbose=verbose);
  verbose && cat(verbose, "Located pathname: ", pathname);

  if (is.null(pathname))
    throw("Failed to located dChip genome information: ", chipType);

  verbose && enter(verbose, "Instantiating ", class(static)[1]);
  res <- newInstance(static, pathname);
  verbose && print(verbose, res);
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
})


setMethodS3("verify", "DChipSnpInformation", function(this, ...) {
  tryCatch({
    df <- readDataFrame(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the dChip SNP information file (",
                                 ex$message, "): ", getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)


setMethodS3("readDataFrame", "DChipSnpInformation", function(this, ...) {
  readFcns <- list(
    "^Mapping10K" = read10K,
    "^Mapping50K" = read50K,
    "^Mapping250K" = read250K
  );

  chipType <- getChipType(this);

  # Try to read with the designated read function.
  res <- NULL;
  for (kk in seq_along(readFcns)) {
    pattern <- names(readFcns)[kk];
    if (regexpr(pattern, chipType) != -1) {
      readFcn <- readFcns[[kk]];
      tryCatch({
        res <- readFcn(this, ...);
      }, error=function(ex) {})
    }
  }

  # If failed, re-try using all read functions.
  if (is.null(res)) {
    for (kk in seq_along(readFcns)) {
      readFcn <- readFcns[[kk]];
      tryCatch({
        res <- readFcn(this, ...);
      }, error=function(ex) {})
      if (!is.null(res))
        break;
    }
  }

  if (is.null(res)) {
    throw("Cannot read dChip SNP information file.  No predefined read function available for this chip type: ", chipType);
  }

  res;
})




setMethodS3("read250K", "DChipSnpInformation", function(this, ..., exclude=c("dbSNP RS ID", "Flank", "FreqAsia", "FreqAfAm", "FreqCauc")) {
  # Example with TABs replaced by semicolons:
  # Probe Set ID;dbSNP RS ID;Flank;Fragment Length Start Stop;FreqAsia;FreqAfAm;FreqCauc
  # SNP_A-1780520;rs16994928;ggatagtgttgacctc[A/G]agtacaggtttcaaaa;496 // 47873735 // 47874230;0.0 ;0.11;0.0
  colClasses <- c(
    "Probe Set ID"="character",
    "dbSNP RS ID"="character",
    "Flank"="character",
    "Fragment Length Start Stop"="character",
    "FreqAsia"="double",
    "FreqAfAm"="double",
    "FreqCauc"="double"
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, private=TRUE)

setMethodS3("read50K", "DChipSnpInformation", function(this, ..., exclude=c("dbSNP RS ID", "Flank", "FreqAsian", "FreqAfAm", "FreqCauc")) {
  colClasses <- c(
    "Probe Set ID"="character",
    "dbSNP RS ID"="character",
    "Flank"="character",
    "Fragment Length Start Stop"="character",
    "FreqAsian"="double",
    "FreqAfAm"="double",
    "FreqCauc"="double"
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, private=TRUE)

setMethodS3("read10K", "DChipSnpInformation", function(this, ..., exclude=c("dbSNP RS ID", "Flank", "Freq Asian", "Freq AfAm", "Freq Cauc"), fill=TRUE) {
  colClasses <- c(
    "Probe Set ID"="character",
    "dbSNP RS ID"="character",
    "Flank"="character",
    "Fragment Length Start Stop"="character",
    "Freq Asian"="double",
    "Freq AfAm"="double",
    "Freq Cauc"="double"
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, fill=fill, ...);
}, private=TRUE)


############################################################################
# HISTORY:
# 2007-01-22
# o Rename argument 'path' to 'rootPath' and added argument 'pattern' to
#   method fromChipType().
# o Just like for genome information file, fromChipType() and readData()
#   were updated to better locate and read dChip SNP information files.
# o Made readData() to support also unknown chip types. That is, if a chip
#   type is not among the hardwired ones, the method will still try to read
#   it using one of the known read functions.  For instance, the genome info
#   file for 10K chips have the same format as the one for the 100K chips.
# 2006-09-17
# o Created from DChipGenomeInformation.R.
############################################################################
