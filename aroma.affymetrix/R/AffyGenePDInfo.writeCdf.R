###########################################################################/**
# @set "class=AffyGenePDInfo"
# @RdocMethod writeCdf
# @alias writeCdf.PDInfoList
# @alias writeCdf.DBPDInfo
#
# @title "Generates an Affymetrix CDF file from a Platform Design (PD) package"
#
# \description{
#   @get "title".
#   Platform Design (PD) packages are also known as "pdInfo" packages.
# }
#
# @synopsis
#
# \arguments{
#  \item{tags}{An optional @character @vector of tags to be added to the CDF
#    filename.}
#  \item{unitsBy}{A @character string specifying how to group units.}
#  \item{path}{The output path where the CDF file is written.
#    If @NULL (default), then it is written to the corresponding
#    \code{annotationData/chipTypes/<chipType>/} directory.}
#  \item{overwrite}{If @TRUE, an existing CDF is overwritten, otherwise
#    an exception is thrown.}
#  \item{verbose}{A @logical or @see "R.utils::Verbose".}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns (invisibly) the pathname to CDF written.
# }
#
# \details{
#   The formal chip type of the CDF is inferred from the AffyGenePDInfo package.
#   The filename of the CDF is constructed from the chip type and any optional
#   tags.
#   To minimize the risk for a corrupt CDF file, the creation of the file
#   is atomic by first writing to a temporary file which is then renamed.
# }
#
# \section{Limitations}{
#   The information available in the PD package is limited and does
#   not contain all information needed to populate a CDF file.
#   In order to workaround these limitations, certain CDF entries
#   are set to predefined/hardwired values.
#   The 'pbase' and 'tbase' entries of the generated CDF file is
#   hardwired to "T" and "A", respectively.
#   Likewise, the 'groupdirection' entry is hardwired to "sense".
# }
#
# \author{
#   Henrik Bengtsson and Guido Hooiveld adopted from \code{pdInfo2Cdf()}
#   written by Samuel Wuest and Mark Robinson.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("writeCdf", "AffyGenePDInfo", function(this, tags=c("*"), unitsBy=c("transcript", "exon"), namesBy=c("fsetid", "id"), path=NULL, overwrite=FALSE, verbose=TRUE, ...) {
  # Early error, iff package is missing
  requireNamespace("DBI") || throw("Package not loaded: DBI")
  dbGetQuery <- DBI::dbGetQuery
  dbListTables <- DBI::dbListTables

  requireNamespace("oligoClasses") || throw("Package not loaded: oligoClasses")
  db <- oligoClasses::db

  requireNamespace("pdInfoBuilder") || throw("Package not loaded: pdInfoBuilder")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  if (is.null(tags)) tags <- c("*");
  tags <- Arguments$getTags(tags, collapse=NULL);

  # Argument 'path':
  if (!is.null(path)) {
    path <- Arguments$getWritablePath(path);
  }

  # Argument 'unitsBy':
  unitsBy <- match.arg(unitsBy);

  # Argument 'namesBy':
  namesBy <- match.arg(namesBy);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # Platform Design object
  pd <- this;


  verbose && enter(verbose, "Generating CDF file from Platform Design (PD) package");
  pkgName <- pd@annotation;
  verbose && cat(verbose, "Platform Design (PD) package: ", pkgName);

  # Infer the chip type
  pkgInfo <- packageDescription(pkgName);
  title <- pkgInfo$Title;
  chipType <- gsub(".* ", "", title);
  # Not needed anymore
  pkgInfo <- title <- NULL;

  # Chip type package name
  pkgNameT <- .cleanPlatformName(chipType);

  # Sanity check
  stopifnot(pkgNameT == pkgName);
  # Not needed anymore
  pkgNameT <- NULL;


  if (is.null(path)) {
    path <- file.path("annotationData", "chipTypes", chipType);
    path <- Arguments$getWritablePath(path);
  }

  verbose && cat(verbose, "Output path: ", path);

  # Asterisk tags?
  if (any(tags == "*")) {
    aTags <- sprintf("by%s-%s", capitalize(unitsBy), namesBy);
    aTags <- c(aTags, pkgName);
    aTags <- Arguments$getTags(aTags);
    tags[tags == "*"] <- aTags;
  }
  tags <- Arguments$getTags(tags);

  chipTypeF <- paste(c(chipType, tags), collapse=",");
  filename <- sprintf("%s.cdf", chipTypeF);
  verbose && cat(verbose, "Filename: ", filename);

  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);
  verbose && cat(verbose, "Pathname to generated CDF: ", pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve chip type dimensions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dim <- .geometry(pd);
  nrows <- dim[1];
  ncols <- dim[2];
  # Not needed anymore
  dim <- NULL;

  nbrOfCells <- nrows*ncols;
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Tags: ", tags);
  verbose && printf(verbose, "Chip type dimensions: %dx%d\n", nrows, ncols);
  verbose && cat(verbose, "Total number of cells (probes): ", nbrOfCells);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving information from PD package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Querying Platform Design database");
  verbose && cat(verbose, "Units by: ", unitsBy);
  verbose && cat(verbose, "Names by: ", namesBy);

  db <- db(pd);
  tables <- dbListTables(db);
  verbose && cat(verbose, "Available database tables:");
  verbose && print(verbose, tables);

  # Sanity checks
  needs <- c("pmfeature", "core_mps");
  if (namesBy == "id") {
    needs <- c(needs, "featureSet");
  }
  for (need in needs) {
    if (!is.element(need, tables)) {
      throw(sprintf("Required table '%s' not in database '%s': %s", need, .annotation(pd), paste(tables, collapse=", ")));
    }
  }

  if (unitsBy == "transcript") {
    # Gives meta-probesets/transcripts
    if (namesBy == "id") {
      sql <- "
        SELECT DISTINCT
--          fid,
--          meta_fsetid,
          featureSet.transcript_cluster_id AS probeset_id,
--          featureSet.exon_id,
          atom,
          x,
          y
        FROM pmfeature
        INNER JOIN core_mps USING(fsetid)
        INNER JOIN featureSet USING(fsetid)
        ORDER BY probeset_id, atom
      ";
    } else if (namesBy == "fsetid") {
      sql <- "
        SELECT DISTINCT
--          fid,
          meta_fsetid AS probeset_id,
          atom,
          x,
          y
        FROM pmfeature
        INNER JOIN core_mps USING(fsetid)
        ORDER BY probeset_id, atom
      ";
    }
  } else if (unitsBy == "exon") {
    # Gives probesets/exons
    if (namesBy == "id") {
      sql <- "
        SELECT DISTINCT
--          fid,
--          fsetid,
--          featureSet.transcript_cluster_id,
          featureSet.exon_id AS probeset_id,
          atom,
          x,
          y
        FROM pmfeature
        INNER JOIN featureSet USING(fsetid)
        ORDER BY probeset_id, atom
      ";
    } else if (namesBy == "fsetid") {
      sql <- "
        SELECT DISTINCT
--          fid,
          fsetid AS probeset_id,
          atom,
          x,
          y
        FROM pmfeature
        ORDER BY probeset_id, atom
      ";
    }
  }

  verbose && cat(verbose, "SQL query:");
  verbose && cat(verbose, sql);

  ff <- dbGetQuery(db, sql);
  verbose && str(verbose, ff);

  # Fix missing values in 'probeset_id'
  id <- ff$probeset_id;
  stopifnot(!is.null(id)); # Sanity check
  nok <- which(is.na(id));
  if (length(nok) > 0L) {
    id <- as.character(id);
    id[nok] <- "<NA>";
    ff$probeset_id <- id;
  }

  # Not needed anymore
  pd <- NULL;  # Not needed anymore
  verbose && str(verbose, ff);
  nbrOfPdCells <- nrow(ff);
  verbose && printf(verbose, "Number of cells (probes) in PD database: %d (%.2f%%) of %d\n", nbrOfPdCells, 100*nbrOfPdCells/nbrOfCells, nbrOfCells);
  verbose && exit(verbose);


  unitNames <- unique(id);
  nbrOfUnits <- length(unitNames);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  verbose && cat(verbose, "Unit names: ", hpaste(unitNames));
  verbose && printf(verbose, "Average number of cells per units: %.2f\n", nbrOfPdCells/nbrOfUnits);

  # Lines speeding up the splitting...
  idGroups <- substr(id, start=1L, stop=4L);
  ffs <- split(ff, f=idGroups);
  names(ffs) <- NULL;
  ffs <- lapply(ffs, FUN=function(u) split(u, f=u$probeset_id));
  ffs <- unlist(ffs, recursive=FALSE);

  attr(ffs, "chipType") <- chipType;
  attr(ffs, "nrows") <- nrows;
  attr(ffs, "ncols") <- ncols;

  class(ffs) <- c("PDInfoList", class(ffs));


  # Sanity checks
  stopifnot(length(ffs) == nbrOfUnits);
  stopifnot(all(is.element(names(ffs), unitNames)));

  verbose && str(verbose, head(ffs, n=3L));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Writing CDF file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- writeCdf(ffs, pathname=pathname, overwrite=overwrite, verbose=less(verbose));

  ## Create checksum
  dfZ <- getChecksumFile(pathname)

  verbose && exit(verbose);

  invisible(pathname);
}) # writeCdf() for AffyGenePDInfo




setMethodS3("writeCdf", "PDInfoList", function(ffs, pathname, overwrite=FALSE, ..., verbose=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pmFeature2List <- function(u) {
    nr <- nrow(u);
    o <- order(u$atom);
    v <- 0L:(nr-1L);
    id <- u$probeset_id[1L];
    g <- list(list(x=u$x[o], y=u$y[o], pbase=rep("T", times=nr), tbase=rep("A", times=nr), atom=v, indexpos=v, groupdirection="sense", natoms=nr, ncellsperatom=1L));
    names(g) <- id;
    # unittype = 'expression' = 1L,
    # cf. help("readCdfUnits", package="affxparser")
    list(unittype=1L, unitdirection=1L, groups=g, natoms=nr, ncells=nr, ncellsperatom=1L, unitnumber=0L);
  } # pmFeature2List()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=!overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);



  verbose && enter(verbose, "Writing PDInfoList to CDF file");
  verbose && cat(verbose, "Pathname to generated CDF: ", pathname);

  nbrOfUnits <- length(ffs);
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);

  chipType <- attr(ffs, "chipType");
  chipType <- Arguments$getCharacter(chipType, length=1L);
  verbose && cat(verbose, "Chip type: ", chipType);

  nrows <- attr(ffs, "nrows");
  nrows <- Arguments$getInteger(nrows, range=c(1L,Inf));
  ncols <- attr(ffs, "ncols");
  ncols <- Arguments$getInteger(ncols, range=c(1L,Inf));
  verbose && cat(verbose, "Number of rows: ", nrows);
  verbose && cat(verbose, "Number of columns: ", ncols);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setting up CDF tree structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Setting up CDF tree structure");

  verbose && enter(verbose, "Setting CDF header");
  cdfHeader <- list(ncols=ncols, nrows=nrows, nunits=nbrOfUnits,
                    nqcunits=0L, refseq="", chiptype=chipType,
                    filename=pathname, rows=nrows,
                    cols=ncols, probesets=nbrOfUnits,
                    qcprobesets=0L, reference="");
  verbose && exit(verbose);

  verbose && enter(verbose, "Setting up CDF units");
  verbose && cat(verbose, "Number of units: ", nbrOfUnits);
  cdfList <- lapply(ffs, FUN=pmFeature2List);
  # Not needed anymore
  ffs <- NULL;  # Not needed anymore
  # Updating 'unitnumber':s
  for (kk in seq_along(cdfList)) {
    cdfList[[kk]]$unitnumber <- kk;
  }
  verbose && str(verbose, head(cdfList, n=3L));

  verbose && cat(verbose, "Unit names: ", hpaste(names(cdfList)));
  verbose && cat(verbose, "Distribution of unit name lengths:");
  verbose && print(verbose, table(nchar(names(cdfList))));

  verbose && exit(verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Writing CDF to file (binary format)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing (binary) CDF file");
  pathname <- cdfHeader$filename;
  verbose && cat(verbose, "Pathname: ", pathname);

  # Write to a temporary file
  if (overwrite && isFile(pathname)) {
    file.remove(pathname);
  }
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

  .writeCdf(pathnameT, cdfheader=cdfHeader, cdf=cdfList,
                  cdfqc=NULL, verbose=verbose, overwrite=overwrite);

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(pathname);
}, protected=TRUE) # writeCdf() for PDInfoList




setMethodS3("writeCdf", "DBPDInfo", function(this, tags=c("*"), unitsBy=c("transcript", "exon"), namesBy=c("fsetid", "id"), path=NULL, overwrite=FALSE, verbose=TRUE, ...) {
  throw(sprintf("writeCdf() for '%s' not implemented: %s", class(this)[1L], .annotation(this)));
})


############################################################################
# HISTORY:
# 2012-12-18 [HB]
# o ROBUSTNESS: Now writeCdf() for PDInfoList sets incremental
#   'unitnumber':s for the CDF units (and no longer as the unit name!).
# o Added argument 'namesBy'.
# 2012-12-17 [HB]
# o Added internal writeCdf() for PDInfoList.
# o GENERALIZATION: Argument 'unitsBy' of writeCdf() for AffyGenePDInfo
#   replaces former argument 'useTranscriptCluster'.
# o Change the default to writeCdf(..., useTranscriptCluster=TRUE) for
#   AffyGenePDInfo.
# o BUG FIX: writeCdf(..., useTranscriptCluster=FALSE) for AffyGenePDInfo
#   could generate empty unit names files for units with 'fsetid' <= 9999.
# 2012-08-16 [MR]
# o Added argument 'useTranscriptCluster' to writeCdf() for AffyGenePDInfo.
# 2012-06-15 [HB]
# o Now writeCdf() for AffyGenePDInfo returns the pathname of the CDF.
# 2011-01-09 [HB]
# o Added writeCdf() for AffyGenePDInfo, which replaces pdInfo2Cdf().
#   An auxillary CEL file is no longer needed to create a CDF from
#   an PDInfo package.  Moreover, contrary pdInfo2Cdf(), the generated
#   CDF now gets a correct/formal Affymetrix chip type.
#
# Below history is for pdInfo2Cdf():
#
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
