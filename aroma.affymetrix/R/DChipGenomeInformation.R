###########################################################################/**
# @RdocClass DChipGenomeInformation
#
# @title "The DChipGenomeInformation class"
#
# \description{
#  @classhierarchy
#
#  This class represents dChip genome information files, which typically
#  contains information about chromosomal locations of the units.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "GenomeInformation".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#   The dChip genome information files for various chip types can be
#   downloaded from \url{http://www.hsph.harvard.edu/cli/complab/dchip/}.
#   Put each file in a
#   directory named identically as the corresponding chip type under the
#   \emph{annotations/} directory, e.g.
#   \emph{annotations/Mapping50K\_Hind240/50k hind genome info AfAm
#   june 05 hg17.xls}.
#   Note that dChip changes the filename and file format slightly between
#   chip types, but currently the @seemethod "byChipType" basically searches
#   for files with names consisting of \code{"genome info"} or
#   \code{"genome_info"}.  At least for the most common chip types, there
#   is no need to rename the files in order for this class to recognize them.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("DChipGenomeInformation", function(...) {
  this <- extend(GenomeInformation(...), "DChipGenomeInformation")
  if (isFile(this)) verify(this)
  this
})


setMethodS3("findByChipType", "DChipGenomeInformation", function(static, chipType, version=NULL, ...) {
  # Argument 'version':
  if (is.null(version))
    version <- ".*";

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pattern <- sprintf("^.*( |_)genome( |_)info(| |_).*%s[.](txt|xls)$",
                                                               version);
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
# @title "Defines a DChipGenomeInformation object by chip type"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{chipType}{A @character string.}
#  \item{version}{An optional @character string specifying the version
#    string, if more than one version is available.}
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "DChipGenomeInformation" object.
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
setMethodS3("byChipType", "DChipGenomeInformation", function(static, chipType, version=NULL, ..., verbose=FALSE) {
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

setMethodS3("verify", "DChipGenomeInformation", function(this, ...) {
  tryCatch({
    df <- readDataFrame(this, nrow=10);
  }, error = function(ex) {
    throw("File format error of the dChip genome information file (",
                                 ex$message, "): ", getPathname(this));
  })
  invisible(TRUE);
}, private=TRUE)


setMethodS3("readDataFrame", "DChipGenomeInformation", function(this, units=NULL, ..., .orderByUnits=FALSE) {
  readFcns <- list(
    "^GenomeWideSNP"  = readGenomeWideSNP,
    "^Mapping10K"     = read50KHg17,
    "^Mapping50K"     = read50KHg17,
    "^Mapping250K"    = read250KHg17,
    "^Mouse430"       = readMouse430,
    "^Mouse430KS"     = readMouse430KenHardwired
  );

  # Get the chip type; this requires that there is a CDF too.
  tryCatch({
    chipType <- getChipType(this);
  }, error = function(ex) {
    print(ex);
  })

  # Validate 'units'?
  if (!is.null(units)) {
    # Locate the CDF
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Try to read with the designated read function.
  res <- NULL;
  for (kk in seq_along(readFcns)) {
    pattern <- names(readFcns)[kk];
    if (regexpr(pattern, chipType) != -1) {
      readFcn <- readFcns[[kk]];
      tryCatch({
        res <- readFcn(this, ...);
      }, error=function(ex) {
        print(ex);
      })
    }
  }

  # If failed, re-try using all read functions.
  if (is.null(res)) {
    for (kk in seq_along(readFcns)) {
      readFcn <- readFcns[[kk]];
      tryCatch({
        res <- readFcn(this, ...);
      }, error = function(ex) {
        print(ex);
      })
      if (!is.null(res))
        break;
    }
  }

  if (is.null(res)) {
    throw("Cannot read dChip annotation data.  No predefined read function available for this chip type: ", chipType);
  }

  # Extract units of interest?
  if (.orderByUnits) {
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    units <- 1:nbrOfUnits(cdf);
  }

  if (!is.null(units)) {
    cdf <- AffymetrixCdfFile$byChipType(chipType);
    unitNames <- getUnitNames(cdf, units=units);
    idxs <- match(unitNames, res[,1]);
    # Not needed anymore
    unitNames <- NULL;
    res <- res[idxs,,drop=FALSE];
    # Not needed anymore
    idxs <- NULL;
  }

  res;
})



setMethodS3("readGenomeWideSNP", "DChipGenomeInformation", function(this, ..., exclude=c("^Strand$", "^dbSNP RS ID$")) {
  colClasses <- c(
    "^Probe Set ID$"="character",
    "^Chromosome$"="character",
    "^Physical Position$"="integer",
    "^Strand$"="character",
    "^dbSNP RS ID$"="character",
    "^$"="NULL"   # Remove empty columns
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, isPatterns=TRUE, na.strings=c("", "---"), ...);
}, private=TRUE)


setMethodS3("read250KHg17", "DChipGenomeInformation", function(this, ..., exclude=c("Expr1002", "Allele A", "dbSNP RS ID")) {
  colClasses <- c(
    "Probe Set ID"="character",
    "Chromosome"="character",
    "Expr1002"="integer",
    "Physical Position"="integer",
    "Allele A"="character",
    "dbSNP RS ID"="character"
  );
  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, private=TRUE)


setMethodS3("read50KHg17", "DChipGenomeInformation", function(this, ..., exclude=c("Genetic Map", "Strand", "Allele A", "dbSNP RS ID", "Freq AfAm", "heterrate")) {
  colClasses <- c(
    "Probe Set ID"="character",
    "Chromosome"="character",
    "Physical Position"="integer",
    "Genetic Map"="",               # Often empty?!
    "Strand"="character",
    "Allele A"="character",
    "dbSNP RS ID"="character",
    "Freq AfAm"="double",
    "heterrate"="double"            # Often empty?!
  );

  readTableInternal(this, pathname=getPathname(this), colClasses=colClasses, exclude=exclude, ...);
}, private=TRUE)


setMethodS3("readMouse430KenHardwired", "DChipGenomeInformation", function(this, ...) {
  colClasses <- c(
    "Probe Set"="character",
    "chromosome"="character",
    "Physical Position"="integer",
    "End"="integer",
    "Strand"="character",
    "Cytoband"="character"
  );

  tableData <- readTableInternal(this, pathname=getPathname(this),
                                                colClasses=colClasses, ...);

  # ad hoc fix: remove "chr" from chromosome
  tableData[,"chromosome"] <- gsub("chr", "", tableData[,"chromosome"]);
  tableData[,"chromosome"] <- gsub("random", "", tableData[,"chromosome"]);

  tableData;
}, private=TRUE)


# @author "KS"
setMethodS3("readMouse430", "DChipGenomeInformation", function(this, ...) {
  colClasses <- c(
    "Probe Set"="character",
    "chromosome"="character",
    "Start"="integer",
    "End"="integer",
    "Strand"="character",
    "Cytoband"="character"
  );

  tableData <- readTableInternal(this, pathname=getPathname(this),
                                                colClasses=colClasses, ...);
  colnames <- colnames(tableData);
  colnames <- gsub("^Start", "Physical Position", colnames);
  colnames(tableData) <- colnames;

  # Remove the "chr" prefix from the "chromosome" elements /KS
  chr <- tableData[,"chromosome"];
  chr <- gsub("chr", "", chr);
  # Ad hoc fix: exclude all chrN_random items /KS
  chr[grep("_random$", chr)] <- NA;
  tableData[,"chromosome"] <- chr;

  tableData;
}, private=TRUE)



############################################################################
# HISTORY:
# 2008-12-29
# o Added argument 'units' to readDataFrame().
# 2007-08-12
# o Added support for dChip's 'snp6.0 genome info hg18.txt' file.  Note,
#   this only contains informations for SNPs, not CN probes.
# 2007-03-23
# o It turns out that KS's readMouse430() was a bit too hardwired.  I've
#   updated it according to his explanations, but I have not tested.
# 2007-03-19 [KS]
# o Added private readMouse430() to the list of read methods that readData()
#   is using to parse DChip genome information files.
# 2007-01-22
# o Rename argument 'path' to 'rootPath' and added argument 'pattern' to
#   method fromChipType().
# o Made fromChipType() identify genome information files more robustly,
#   e.g. the exact chip type does not have to be part of the prefix.
#   Indeed, it now accepts the default dChip filenames for the 10K, the
#   50K and the 250K SNP chips.
# o Made readData() to support also unknown chip types. That is, if a chip
#   type is not among the hardwired ones, the method will still try to read
#   it using one of the known read functions.  For instance, the genome info
#   file for 10K chips have the same format as the one for the 100K chips.
# 2006-09-15
# o Created from DChip.R in old(!) aroma.snp.
# 2005-11-15
# o Added support for skipping header in readSampleInformationFile().
# 2005-10-31
# o Created.
############################################################################
