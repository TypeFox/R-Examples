###########################################################################/**
# @set "class=AffymetrixCdfFile"
# @RdocMethod createExonByTranscriptCdf
#
# @title "Creates an exon-by-transcript CDF"
#
# \description{
#  @get "title" based on the probesets defined in an "exon-only" CDF
#  and transcript-exon mapping of a NetAffx probeset annotation data file.
# }
#
# @synopsis
#
# \arguments{
#  \item{cdf}{An @see "aroma.affymetrix::AffymetrixCdfFile" specifying
#     an "exon-only" CDF, which defines the exon-specific probesets
#     that will go into the new CDF. For more details, see below.}
#  \item{csv}{An @see "aroma.affymetrix::AffymetrixNetAffxCsvFile"
#     specifying the Affymetrix NetAffx CSV probeset annotation file
#     that contains the transcript-exon mapping.}
#  \item{tags}{Additional tags added to the filename of created CDF,
#      i.e. <chiptype>,<tags>.cdf.}
#  \item{path}{The output path where the custom CDF will be written.}
#  \item{type}{A @character string specifying the type of CDF to be written.}
#  \item{subsetBy}{An optional @character specifying the name of a column
#     in the annotation file to subset against.  The column will be parsed
#     as the data type of argument \code{within}.}
#  \item{within}{A @vector of values accepted for the \code{subsetBy} column.}
#  \item{...}{Additional arguments passed to \code{readDataFrame()} of
#     @see "aroma.affymetrix::AffymetrixNetAffxCsvFile", e.g. \code{nrow}.}
#  \item{overwrite}{If @TRUE, an existing CDF is overwritten.}
#  \item{verbose}{...}
# }
#
# \value{
#   Returns an @see "aroma.affymetrix::AffymetrixCdfFile".
# }
#
# \section{Requirements for the "exon-only" CDF}{
#   The template CDF - argument \code{cdf} - should be an "exon-only" CDF:
#   each unit has one group/probeset, which is the exon.
#   An example of this is the "unsupported" HuEx-1_0-st-v2.cdf
#   from Affymetrix, which has 1,432,154 units.
#   Such "exon-only" CDFs do not contain information about clustering
#   exons/probesets into gene transcripts.
#   The CDF may also contain a number of non-exon probesets corresponding
#   to control probes, which can contain \emph{very} large numbers of
#   probes per probeset. Such units are dropped/ignored by this method.
# }
#
# \section{Ordering of transcripts and exons within transcripts}{
#   The transcripts (=units) will be ordered as they appear in the
#   NetAffx annotation file.
#   Within each transcript (=unit), the exons (=groups) are ordered
#   lexicographically by their names.
#   %% (Before Jan 28, 2008 (rev 3911) sorting was not done).
# }
#
# \section{Naming of transcripts and exons}{
#   In the created CDF, each unit corresponds to one transcript cluster,
#   and each group within a unit corresponds to the exons within
#   that transcript cluster.  Thus, the unit names correspond to the
#   transcript cluster names and the group names correspond to the
#   exon names.
#
#   The exon names are defined by unit names of the exon-only CDF,
#   whereas the transcript names are defined by the
#   \code{transcriptClusterId} column in the NetAffx annotation data file.
#   These transcript and exon names are often "non-sense" integers.
#   In order to map these to more genome-friendly names, use the various
#   annotations available in the NetAffx annotation data file.
# }
#
# \examples{\dontrun{
# # The exon-only CDF
# cdf <- AffymetrixCdfFile$byChipType("HuEx-1_0-st-v2");
#
# # The NetAffx probeset annotation data file
# csv <- AffymetrixNetAffxCsvFile("HuEx-1_0-st-v2.na24.hg18.probeset.csv", path=getPath(cdf));
#
# # Create a CDF containing all core probesets:
# cdfT <- createExonByTranscriptCdf(cdf, csv=csv, tags=c("*,HB20110911"));
# print(cdfT);
#
# # Create CDF containing the core probesets with 3 or 4 probes:
# cdfT2 <- createExonByTranscriptCdf(cdf, csv=csv,
#             tags=c("*,bySize=3-4,HB20110911"),
#             subsetBy="probeCount", within=c("3", "4"));
# print(cdfT2);
# }}
#
# \author{
#   Henrik Bengtsson adopted from \code{createTranscriptCDF()} written
#   by Ken Simpson, Elizabeth Purdom and Mark Robinson.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("createExonByTranscriptCdf", "AffymetrixCdfFile", function(cdf, csv, tags=c("*"), path=getPath(cdf), type=c("all", "core", "extended", "full", "main", "control", "cds"), subsetBy=NULL, within=NULL, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  inMemory <- FALSE; # Turns out not to make a big difference. /HB 2011-09-10
  if (inMemory) {
    cdfTree <- NULL;
    getCdfUnits <- function(units) {
      if (is.null(cdfTree)) {
        verbose && enter(verbose, "Reading the complete CDF");
        cdfTree <<- .readCdf(cdfPathname);
        verbose && exit(verbose);
      }
      cdfTree[units];
    } # getCdfUnits()
  } else {
    getCdfUnits <- function(units) {
      ## SLOW iff units by units!  /HB 2011-09-09
      # (Actually, not that much slower. /HB 2011-09-10)
      .readCdf(cdfPathname, units=units);
    } # getCdfUnits()
  }

  makeTranscriptCdfUnit <- function(transcriptName) {
    # Identify all probesets (=exons) part of the transcript of interest
    keep <- which(psData[,"transcriptClusterId"] == transcriptName);
    groupNames <- psData[keep,"probesetId"];
    groupNames <- sort(groupNames);

    # Find those exons (=units) in the "exon-only" CDF
    groupUnits <- indexOf(cdf, names=groupNames);  ## <= FAST - cached!
    # Sanity check
    stopifnot(all(is.finite(groupUnits)));

    # Read the CDF units of those exons
    cdfList <- getCdfUnits(units=groupUnits);

    # Extract the groups
    groups <- lapply(cdfList, FUN=.subset2, "groups");
    # Sanity check (of the assumption of single group exon units)
    nGroups <- sapply(groups, FUN=length);
    stopifnot(all(nGroups == 1));
    groups <- lapply(groups, FUN=.subset2, 1);

    # Identify the unit type
    unittypes <- sapply(cdfList, FUN=.subset2, "unittype");
    unittype <- unittypes[1];
    # Sanity check (of assumption)
    stopifnot(all(unittypes == unittype));

    # Identify the unit direction
    unitdirections <- sapply(cdfList, FUN=.subset2, "unitdirection");
    unitdirection <- unitdirections[1];
    # Sanity check (of assumption)
    stopifnot(all(unitdirections == unitdirection));

    # Identify the number of cells per atom
    ncellsperatoms <- sapply(cdfList, FUN=.subset2, "ncellsperatom");
    ncellsperatom <- ncellsperatoms[1];
    # Sanity check (of assumption)
    stopifnot(all(ncellsperatoms == ncellsperatom));

    # Identify the number of atoms
    natoms <- lapply(cdfList, FUN=.subset2, "natoms");
    natoms <- sum(unlist(natoms, use.names=FALSE));

    # Identify the number of cells
    ncells <- lapply(cdfList, FUN=.subset2, "ncells");
    ncells <- sum(unlist(ncells, use.names=FALSE));

    # Return a CDF unit
    list(
      groups=groups,
      unittype=unittype,
      unitdirection=unitdirection,
      natoms=natoms,
      ncells=ncells,
      ncellsperatom=ncellsperatom,
      unitnumber=unitnumber
    );
  } # makeTranscriptCdfUnit()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'type':
  type <- match.arg(type);

  nbrOfCdfUnits <- nbrOfUnits(cdf);

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- unlist(strsplit(tags, split=",", fixed=TRUE));
    tags <- tags[nzchar(tags)];
  }

  # Argument 'csv':
  if (!inherits(csv, "AffymetrixNetAffxCsvFile")) {
    pathname <- csv;
    csv <- AffymetrixNetAffxCsvFile(pathname);
  }

  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Creating exon-by-transcript CDF");

  cdfPathname <- getPathname(cdf);
  chipType <- getChipType(cdf);
  verbose && cat(verbose, "\"Exon-only\" CDF: ", cdfPathname);
  verbose && cat(verbose, "Chip type: ", chipType);
  verbose && cat(verbose, "Number of units in template CDF: ", nbrOfCdfUnits);

  # Insert asterisk tags
  if (any(tags == "*")) {
    idxs <- which(tags == "*");
    genomeTag <- getGenomeBuild(csv);
    naTagA <- sprintf("na%s", getNetAffxBuild(csv));
##    naTagB <- format(getNetAffxDate(csv), format="%Y%m%d");
##    naTagB <- sprintf("na%s", naTagB);
    asteriskTags <- c(type, naTagA, genomeTag);
    asteriskTags <- asteriskTags[nzchar(asteriskTags)];
    asteriskTags <- paste(asteriskTags, collapse=",");
    verbose && cat(verbose, "Asterisk tags: ", asteriskTags);
    tags[idxs] <- asteriskTags;
  }

  # Validate the output pathname already here
  chipTypeF <- paste(c(chipType, tags), collapse=",");
  verbose && cat(verbose, "Tags: ", tags);
  verbose && cat(verbose, "Chip type (fullname): ", chipTypeF);
  filename <- sprintf("%s.cdf", chipTypeF);
  pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=!overwrite);
  # Not needed anymore
  filename <- NULL;

  verbose && cat(verbose, "NetAffx annotation data file");
  verbose && print(verbose, csv);

  # set up the comparisons to do for paring down
  if (is.null(within)) {
    within <- switch(type,
      "core"     = "core",
      "extended" = c("core", "extended"),
      "full"     = c("core", "extended", "full"),
      "main"     = "main",
      "control"  = c("control->affx", "control->chip",
                     "control->bgp->antigenomic", "control->bgp->genomic",
                     "normgene->exon", "normgene->intron"),
      "all"      = NULL,
      "cds"      = "1"
    );
  }


  if (is.null(subsetBy)) {
    if (type %in% c("core", "extended", "full")) {
      subsetBy <- "level";
    }
    if (type %in% c("main", "control")) {
      subsetBy <- "probesetType";
    }
    if (type %in% "cds") {
      subsetBy <- "hasCds";
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read NetAffx CSV file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading the NetAffx annotation file");
  verbose && print(verbose, csv);

  subsetByPattern <- storage.mode(within);
  names(subsetByPattern) <- subsetBy;
  colClasses <- c("(probesetId|transcriptClusterId)"="character", subsetByPattern);
  verbose && cat(verbose, "Column class patterns:");
  verbose && print(verbose, colClasses);

  # Read CSV data
  psData <- readDataFrame(csv, colClasses=colClasses, ...);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify subset of probesets to include
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying subset of probesets to include");
  nbrOfProbesets <- nrow(psData);

  verbose && cat(verbose, "Number of CSV probesets: ", nbrOfProbesets);

  if (nbrOfProbesets > nbrOfCdfUnits) {
    warning(sprintf("Annotation has %d lines (probesets); original CDF only has %d.", nbrOfProbesets, nbrOfCdfUnits));
  }

  verbose && enter(verbose, "Keeping probesets with names matching the CDF unit names");
  verbose && cat(verbose, "Number of CDF units: ", nbrOfCdfUnits);

  # Drop probesets in annotation data file that does not exist in the CDF
  units <- indexOf(cdf, names=psData[,"probesetId"]);
  if (any(is.na(units))) {
    warning(sprintf("%d probesets in annotation are not in original CDF; will be discarded.", length(which(is.na(units)))));
    psData <- psData[!is.na(units),];
    nbrOfProbesets <- nrow(psData);
  }
  # Not needed anymore
  units <- NULL;
  verbose && cat(verbose, "Number of (filtered) CSV probesets: ", nbrOfProbesets);
  verbose && exit(verbose);


  if (nbrOfProbesets < nbrOfCdfUnits) {
    verbose && printf(verbose, "Fewer (filtered) probesets in annotation than orginal CDF (%d compared to %d)\n", nbrOfProbesets, nbrOfCdfUnits);
  }

  if (!is.null(within)) {
    verbose && enter(verbose, "Subsetting probesets");
    verbose && printf(verbose, "Keeping probesets whose '%s' is in (%s)\n", subsetBy, hpaste(within, maxHead=Inf));
    keep <- which(is.element(psData[,subsetBy], within));
    psData <- psData[keep,];
    verbose && printf(verbose, "%d probesets match requirement (out of %d valid probesets in annotation)\n", length(keep), nbrOfProbesets);
    nbrOfProbesets <- nrow(psData);

    # Nothing todo?
    if (nbrOfProbesets == 0) {
      throw("Cannot create CDF. No probesets remaining.");
    }

    # Not needed anymore
    keep <- NULL;
    verbose && exit(verbose);
  }
  # Not needed anymore
  nbrOfProbesets <- NULL;
  gc();

  # Sanity check
  units <- indexOf(cdf, names=psData[,"probesetId"]);
  stopifnot(all(is.finite(units)));
  # Not needed anymore
  units <- NULL;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build custom CDF list structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Grouping exons into transcripts");

  transcriptNames <- unique(psData[,"transcriptClusterId"]);
  nbrOfTranscripts <- length(transcriptNames);
  nbrOfExons <- nrow(psData);

  verbose && cat(verbose, "Number of transcripts: ", nbrOfTranscripts);
  verbose && cat(verbose, "Number of exons: ", nbrOfExons);
  verbose && printf(verbose, "Average number of exons/transcript: %.2f\n", nbrOfExons/nbrOfTranscripts);

  if (nbrOfTranscripts > 100) {
    verbose && printf(verbose, "Number of transcripts done so far:\n");
  }

  # Allocate the CDF list structure
  cdfList <- vector("list", length=nbrOfTranscripts);
  names(cdfList) <- transcriptNames;

  unitnumber <- 0L;
  for (jj in seq_len(nbrOfTranscripts)) {
    unitnumber <- unitnumber + 1L;
    if (unitnumber %% 100 == 0) {
      if (isVisible(verbose)) printf(", %d", unitnumber);
    }
    cdfList[[jj]] <- makeTranscriptCdfUnit(transcriptNames[jj]);
  } # for (jj ...)
  verbose && printf(verbose, "\n");

  # Sanity check
  nbrOfExons2 <- sum(sapply(cdfList, FUN=function(unit) length(unit$group)));
#  verbose && cat(verbose, "Number of exons included: ", nbrOfExons);
  stopifnot(nbrOfExons2 == nbrOfExons);

  # Not needed anymore
  # Not needed anymore
  psData <- NULL;

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write custom CDF to file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing CDF");

  verbose && cat(verbose, "Chip type: ", chipTypeF);
  verbose && cat(verbose, "Number of units (\"transcripts\"): ", nbrOfTranscripts);

  verbose && enter(verbose, "Setting up CDF header");
  # Copy header and QC info from original CDF
  qc <- .readCdfQc(cdfPathname);
  hdr <- .readCdfHeader(cdfPathname);
  hdr$chiptype <- chipTypeF;
  hdr$filename <- pathname;
  hdr$probesets <- nbrOfExons;
  hdr$nunits <- nbrOfTranscripts;

  verbose && exit(verbose);


  verbose && enter(verbose, "Writing to file");

  # Remove existing file?
  if (overwrite && isFile(pathname)) {
    file.remove(pathname);
  }

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);
  .writeCdf(pathnameT, cdfheader=hdr, cdf=cdfList, cdfqc=qc, overwrite=TRUE, verbose=10);

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);
  verbose && exit(verbose);

  verbose && exit(verbose);  # "Writing CDF...exit"

  cdfT <- newInstance(cdf, pathname);
  verbose && cat(verbose, "Created transcript-by-exon CDF:");
  verbose && print(verbose, cdfT);

  verbose && exit(verbose);

  cdfT;
}) # createExonByTranscriptCdf()


############################################################################
# HISTORY:
# 2011-09-11 [HB]
# o DOCUMENTATION: Improved the help on createExonByTranscriptCdf().
# o Renamed to createExonByTranscriptCdf() (from createTranscriptCDF()).
# o ROBUSTNESS: Now createTranscriptCDF() uses the AffymetrixNetAffxCsvFile
#   class to read the annotation data.
# o CLEANUP: User no longer have to specify the indices of the annotation
#   data columns to be read. Addressing is instead done by column names.
# o CLEANUP: Dropped argument 'writeCdf'.
# o CLEANUP: Dropped also argument 'nrow'.  Now the list of arguments
#   are rather clean.
# 2011-09-10 [HB]
# o CLEANUP: Now createTranscriptCDF() no longer adds a dummy QC units,
#   which used to be a workaround for a bug in affxparser that was
#   corrected in affxparser v1.7.0 (2006-10-25).
# o Now argument 'type' lists all possible values.
# o Renamed internal getCdfInfo() to makeTranscriptCdfUnit().
# o ROBUSTNESS: The internal getCdfInfo() assumed that the exon CDF units
#   have all the same (i) unit types, (ii) unit directions, and (iii)
#   number of cells per atom.  Now, if CDF that does not meet those
#   assumption is used, an exception is thrown.
# 2011-09-09 [HB]
# o SPEED UP: Added argument 'inMemory' to allow to read the CDF structure
#   once and cache it in memory, instead of read transcript by transcript.
# o DOCUMENTATION: Created a good stubb of an Rdoc help.
# o ROBUSTNESS: Now the CDF is written atomically (via a temporary file).
# o NOTE: The createTranscriptCDF() sorts the groups (=exons)
#   lexicographically by name.  This was note the case when
#   'HuEx-1_0-st-v2,R3,A20071112,EP.zip' was generated (in December 2007).
# o SPEED UP: Using unlist(...,, use.names=FALSE).
# o Replaced all cat() with Verbose statements.
# o RCC: Renamed all iterators to double-letters, e.g. 'jj' instead of 'j'.
# o CLEANUP: Replaced all cat(paste()) with printf().
# o CLEANUP: Replaced all stop(paste()) with throw().
# o ROBUSTNESS: Now arguments are validated using the Arguments class.
# 2011-09-01 [HB]
# o Now createTranscriptCDF() is a method for the AffymetrixCdfFile class.
# o Tidying up the code for readability, e.g. adding spaces, semicolons,
#   adjusting indentation etc.
############################################################################
