###########################################################################/**
# @RdocFunction pdInfo2Cdf
#
# @title "Generates an Affymetrix CDF file from a Platform Design (PD) package and a auxillary CEL file for the same chip type"
#
# \description{
#   @get "title".
#   Platform Design (PD) packages are also known as "pdInfo" packages.
#
#   \emph{Disclaimer: This is a user-contributed function.}
#
#   \emph{Instead of using this method, we recommend to use
#   \code{\link[=writeCdf.AffyGenePDInfo]{writeCdf}()}
#   for the \code{AffyGenePDInfo} class.}
# }
#
# @synopsis
#
# \arguments{
#  \item{pdpkg}{A @character string for an existing PD package.}
#  \item{celfile}{The pathname to an auxillary CEL for the same chip type.}
#  \item{overwrite}{If @TRUE, an existing CDF is overwritten, otherwise
#    an exception is thrown.}
#  \item{verbose}{A @logical or @see "R.utils::Verbose".}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the pathname to CDF written.
#   The CDF filename is generated from the name of the PD package.
# }
#
# \section{Limitations}{
#   The information available in the PD package is limited and does
#   not contain all information needed to populate a CDF file.
#   In order to workaround these limitations, certain CDF entries
#   are set to predefined/hardwired values.
#   The 'pbase' and 'tbase' entries of the generated CDF file is
#   hardwired to "T" and "A", respectively.  Likewise, the 'groupdirection'
#   entry is hardwired to "sense".
# }
#
# \author{
#   Maintained by Mark Robinson.
#   Original code by Samuel Wuest.
#   Code improvements by Henrik Bengtsson.
# }
#
# \seealso{
#   Instead of using this method, we recommend to use
#   \code{\link[=writeCdf.AffyGenePDInfo]{writeCdf}()}
#   for the \code{AffyGenePDInfo} class.
# }
#
# @keyword internal
#*/###########################################################################
pdInfo2Cdf <- function(pdpkg, celfile, overwrite=FALSE, verbose=TRUE, ...) {
  requireNamespace("oligo") || throw("Package not loaded: oligo")
  read.celfiles <- oligo::read.celfiles

  requireNamespace("DBI") || throw("Package not loaded: DBI")
  dbGetQuery <- DBI::dbGetQuery

  requireNamespace("oligoClasses") || throw("Package not loaded: oligoClasses")
  db <- oligoClasses::db

  .require <- require
  .require("pdInfoBuilder") || throw("Package not loaded: pdInfoBuilder")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pmFeature2List <- function(u) {
    nr <- nrow(u);
    o <- order(u$atom);
    v <- 0:(nr-1);
    id <- u$fsetid[1];
    g <- list(list(x=u$x[o], y=u$y[o], pbase=rep("T", nr), tbase=rep("A", nr),
              atom=v, indexpos=v, groupdirection="sense", natoms=nr,
              ncellsperatom=1));
    names(g) <- id;
    list(unittype=1, unitdirection=1, groups=g, natoms=nr, ncells=nr,
         ncellsperatom=1, unitnumber=id);
  } # pmFeature2List()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pdpkg':
  pdpkg <- Arguments$getCharacter(pdpkg);

  # Argument 'celfile':
  celfile <- Arguments$getReadablePathname(celfile, mustExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);


  verbose && enter(verbose, "Generating CDF file from Platform Design (PD) package");
  verbose && cat(verbose, "Platform Design (PD) package: ", pdpkg);

  pdName <- gsub("\\.", "", pdpkg);

  filename <- sprintf("%s.cdf", pdName);
  filename <- Arguments$getWritablePathname(filename, mustNotExist=!overwrite);

  verbose && cat(verbose, "CDF file to be generated: ", filename);

  # Loading the required PD package.
  .require(pdpkg, character.only=TRUE) ||
                 throw("Platform Design (PD) package not loaded: ", pdpkg);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving information from the CEL file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading auxillary CEL file");
  verbose && cat(verbose, "Pathname: ", celfile);

  verbose && enter(verbose, "Reading CEL file header");
  hdr <- .readCelHeader(celfile);
  nrows <- as.integer(hdr$rows);
  ncols <- as.integer(hdr$cols);
  chipType <- hdr$chiptype;
  # Not needed anymore
  hdr <- NULL;  # Not needed anymore
  nbrOfCells <- nrows*ncols;
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && printf(verbose, "Chip type dimensions: %dx%d\n", nrows, ncols);
  verbose && cat(verbose, "Total number of cells (probes): ", nbrOfCells);
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading complete CEL file");
  cel <- read.celfiles(filenames=celfile, pkgname=pdpkg);
  verbose && exit(verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving information from PD package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving Platform Design database");
  pd <- .getPlatformDesign(cel);
  # Not needed anymore
  cel <- NULL;  # Not needed anymore
  verbose && exit(verbose);

  verbose && enter(verbose, "Querying Platform Design database");
  ff <- dbGetQuery(db(pd), "select * from pmfeature");
  # Not needed anymore
  pd <- NULL;  # Not needed anymore
  verbose && str(verbose, ff);
  nbrOfPdCells <- nrow(ff);
  verbose && printf(verbose, "Number of cells (probes) in PD database: %d (%.2f%%) of %d\n",
                    nbrOfPdCells, 100*nbrOfPdCells/nbrOfCells, nbrOfCells);
  verbose && exit(verbose);

  verbose && enter(verbose, "Creating list from query table");
  # three 3 lines speed up the splitting ...
  ffs <- split(ff, substr(ff$fsetid, start=1, stop=4));
  ffs <- lapply(ffs, FUN=function(u) split(u, u$fsetid));
  ffs <- unlist(ffs, recursive=FALSE);
  names(ffs) <- substr(names(ffs), start=6, stop=nchar(names(ffs)));
  nbrOfUnits <- length(ffs);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  verbose && printf(verbose, "Average number of cells per units: %.2f\n", nbrOfPdCells/nbrOfUnits);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up CDF tree structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up CDF tree structure");

  verbose && enter(verbose, "Setting CDF header");
  newCdfHeader <- list(ncols=ncols, nrows=nrows, nunits=nbrOfUnits,
                       nqcunits=0, refseq="", chiptype=pdName,
                       filename=filename, rows=nrows,
                       cols=ncols, probesets=nbrOfUnits,
                       qcprobesets=0, reference="");
  verbose && exit(verbose);

  verbose && enter(verbose, "Setting up CDF units");
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  newCdfList <- lapply(ffs, FUN=pmFeature2List);
  # Not needed anymore
  ffs <- NULL;  # Not needed anymore
  verbose && exit(verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Writing CDF to file (binary format)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing (binary) CDF file");
  pathname <- newCdfHeader$filename;
  verbose && cat(verbose, "Pathname: ", pathname);
  res <- .writeCdf(pathname, cdfheader=newCdfHeader, cdf=newCdfList,
                  cdfqc=NULL, verbose=verbose, overwrite=overwrite);
  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
} # pdInfo2Cdf()



############################################################################
# HISTORY:
# 2012-03-23
# o Now pdInfo2Cdf() helps 'R CMD check' to locate read.celfiles()
#   by explicitly requiring 'oligo'.
# 2010-12-04 [HB]
# o Added more verbose output.
# o DOCUMENTATION: Added more Rd documentation.
# o BUG FIX: Local variable 'pdName' of pdInfo2Cdf() was used before it
#   was defined.  Thanks to Guido Hooiveld at the Wageningen University,
#   Netherlands, for reporting this.
# 2010-05-20 [HB]
# o Renamed PdInfo2Cdf() to pdInfo2Cdf().  Keeping old one for backward
#   compatibility for a while.
# 2010-05-19 [HB]
# o BUG FIX: PdInfo2Cdf() would write dimension (rows,rows) in the CDF
#   header instead of (rows,cols).  Thanks Kasper Daniel Hansen for
#   reporting this.
# 2009-10-16 [HB]
# o MEMORY OPTIMIZATION: Cleaning out non needed objects sooner.
# o Added verbose a'la R.utils.
# o Added some validation of arguments.
# o Tidied up the code structure.
# 2009-01-13 [MR]
# o Added. "This script has been written to generate a .cdf-file from an
#   "pd.XXXX" package, such as those build with pdInfoBuilder.
#   The original was written by Samuel Wuest, modified by Mark Robinson
#   (around 12 Jan 2009) to be generic."
# 2008-??-?? [SW]
# o Created by Samuel Wuest.
############################################################################
