###########################################################################/**
# @RdocDefault bpmapCluster2Cdf
#
# @title "Creates a CDF from tiling-array BPMAP file"
#
# \description{
#   @get "title".\cr
#
#   \emph{
#    NOTE: This method applies only to Affymetrix tiling arrays!
#    Furthermore, it is likely to be more useful for promoter tiling arrays
#    and less so for whole-genome tiling arrays.
#   }
# }
#
# @synopsis
#
# \arguments{
#  \item{pathname}{The pathname to the BPMAP file.}
#  \item{chipType, tags}{The chip type and optional tags of the CDF to
#    be written.}
#  \item{rows, cols}{Two positive @integers specifying the probe dimension
#    of the chip.  It is important to get this correct.  They can be
#    inferred from the CEL header of a CEL file for this chip,
#    cf. @see "affxparser::readCelHeader".}
#  \item{maxProbeDistance}{A positive @integer specifying the maximum
#    genomic distance (in basepairs) allowed between two probes in order
#    to "cluster" those two probes into the same CDF units.  Whenever the
#    distance is greater, the two probes end up in two different CDF units.}
#  \item{minNbrOfProbes}{A positive @integer specifying the minimum number
#    of probes required in a CDF unit.  If fewer, those probes are dropped.}
#  \item{groupName}{A @character string specifying which BPMAP sequences
#     to keep.  Sequence with this group name is kept, all others are
#     excluded.}
#  \item{field}{A @character string.}
#  \item{stringRemove}{An (optional) regular expression.}
#  \item{...}{Optional arguments passed to @see "affxparser::readBpmap".}
#  \item{flavor}{Specifying which version of BPMAP-to-CDF generator
#     to use. The default is always to use the most recent one, which
#     is also the recommended one.  Previous versions are kept only for
#     backward compatibility (and may be dropped at anytime).}
#  \item{path}{The directory where the CDF file will be written.
#     If \code{"*"} (default), it will be written to
#     \code{annotationData/chipTypes/<chipType>/}.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns (invisibly) a the pathname of the created CDF file.
#  The created CDF is written to the current directory.
# }
#
# \details{
#   This method applies only to Affymetrix tiling arrays.  It is likely
#   to be useful for promoter tiling arrays and less so for whole-genome
#   tiling arrays.
#   Flavor \code{"v2"} replaced \code{"v1"} as aroma.affymetrix v2.5.4
#   (June 21, 2012). For details, see \code{news(package="aroma.affymetrix")}.
# }
#
# \author{
#   Henrik Bengtsson adopted from Mark Robinson standalone/online version
#   as of July 11, 2011.
# }
#
# @keyword "internal"
#*/###########################################################################
setMethodS3("bpmapCluster2Cdf", "default", function(pathname, chipType, tags=NULL, rows, cols, maxProbeDistance=3000L, minNbrOfProbes=30L, groupName=gsub("_.*", "", chipType), field="fullname", stringRemove=sprintf("%s:.*;", groupName), ..., flavor=c("v2", "v1"), path="*", verbose=-10) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getRangeX <- function(bpmapList) {
    pos <- lapply(bpmapList, FUN=.subset, c("mmx", "pmx"));
    pos <- unlist(pos, use.names=FALSE);
    range(pos, na.rm=TRUE);
  } # getRangeX()

  getRangeY <- function(bpmapList) {
    pos <- lapply(bpmapList, FUN=.subset, c("mmy", "pmy"));
    pos <- unlist(pos, use.names=FALSE);
    range(pos, na.rm=TRUE);
  } # getRangeY()

  getGroupNames <- function(bpmapList) {
    names <- sapply(bpmapList, FUN=function(seq) {
      seqInfo <- seq$seqInfo;
      seqInfo$groupname;
    });
    # Sanity check
    stopifnot(length(names) == length(bpmapList));
    names;
  } # getGroupNames()

  getNumberOfProbes <- function(bpmapList) {
    sapply(bpmapList, FUN=function(u) u$seqInfo$numberOfHits);
  } # getNumberOfProbes()

  getStartPositions <- function(bpmapList) {
    lapply(bpmapList, FUN=function(u) u$startpos);
  } # getStartPositions()

  bpmapUnit2df <- function(u) {
    mmx <- u[["mmx"]];
    mmy <- u[["mmy"]];
    mmx <- if (all(mmx == 0L) || is.null(mmx)) 0L else mmx;
    mmy <- if (all(mmy == 0L) || is.null(mmy)) 0L else mmy;

    seqInfo <- u$seqInfo;
    res <- data.frame(seqname=seqInfo[[field]], groupname=seqInfo$groupname, u[c("pmx", "pmy")], mmx=mmx, mmy=mmy, u[c("probeseq", "strand", "startpos", "matchscore")], stringsAsFactors=FALSE);

    # Order by start positions
    o <- order(u[["startpos"]]);
    res <- res[o,,drop=FALSE];

    res;
  } # bpmapUnit2df()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validating arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  pathname <- Arguments$getReadablePathname(pathname);

  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
  }

  # Argument 'minNbrOfProbes':
  minNbrOfProbes <- Arguments$getInteger(minNbrOfProbes, range=c(1,Inf));

  # Argument 'maxProbeDistance':
  maxProbeDistance <- Arguments$getInteger(maxProbeDistance, range=c(1,Inf));

  # Argument 'rows':
  rows <- Arguments$getInteger(rows, range=c(1,Inf));

  # Argument 'cols':
  cols <- Arguments$getInteger(cols, range=c(1,Inf));

  # Argument 'groupName':
  groupName <- Arguments$getCharacter(groupName);

  # Argument 'field':
  field <- Arguments$getCharacter(field);

  # Argument 'stringRemove':
  if (!is.null(stringRemove)) {
    stringRemove <- Arguments$getCharacter(stringRemove);
  }

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'path':
  if (is.null(path)) {
    path <- ".";
  } else {
    path <- Arguments$getCharacter(path);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);



  verbose && enter(verbose, "Generating CDF from BPMAP");
  chipTypeF <- paste(c(chipType, tags), collapse=",");
  verbose && cat(verbose, "Full chip type: ", chipTypeF);
  chipType <- gsub(",.*", "", chipTypeF);

  if (path == "*") {
    path <- file.path("annotationData", "chipTypes", chipType);
  }
  path <- Arguments$getWritablePath(path);
  verbose && cat(verbose, "Output path: ", path);


  cdfFilename <- sprintf("%s.cdf", chipTypeF);
  cdfFilename <- Arguments$getFilename(cdfFilename);
  verbose && cat(verbose, "CDF pathname: ", cdfFilename);
  cdfPathname <- Arguments$getWritablePathname(cdfFilename, path=path, mustNotExist=TRUE);

  ppsPathname <- sprintf("%s.pps", chipTypeF);
  ppsPathname <- Arguments$getFilename(ppsPathname);
  verbose && cat(verbose, "PPS pathname: ", ppsPathname);
  ppsPathname <- Arguments$getWritablePathname(ppsPathname, path=path, mustNotExist=TRUE);

  verbose && enter(verbose, "Reading BPMAP file");
  verbose && cat(verbose, "Pathname: ", pathname);
  hdr <- .readBpmapHeader(pathname);
  verbose && cat(verbose, "File version: ", hdr$version);
  verbose && cat(verbose, "Number of sequences: ", hdr$numSequences);
  bpmapList <- .readBpmap(pathname, readMatchScore=TRUE, ...);
  nbrOfSeqs <- length(bpmapList);
  # Sanity check
  stopifnot(nbrOfSeqs <= hdr$numSequences);
  # Not needed anymore
  hdr <- NULL; # Not needed anymore
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate the chip dimensions further.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  maxY <- getRangeY(bpmapList)[2] + 1L;
  if (maxY > rows) {
    throw("Argument 'rows' is too small. There exist probes with a Y position that is greater: ", maxY, " > ", rows);
  }

  maxX <- getRangeX(bpmapList)[2] + 1L;
  if (maxX > cols) {
    throw("Argument 'cols' is too small. There exist probes with an X position that is greater: ", maxX, " > ", cols);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Filtering sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying sequences to be used in CDF");

  verbose && enter(verbose, "Extracting sequences with group names of interest");
  verbose && cat(verbose, "Group names to keep: ", groupName);
  nbrOfSeqs0 <- nbrOfSeqs;
  verbose && cat(verbose, "Number of sequences before: ", nbrOfSeqs);
  verbose && cat(verbose, "Number of probes before: ", sum(getNumberOfProbes(bpmapList)));

  groupNames <- getGroupNames(bpmapList);
  verbose && cat(verbose, "Group names: ", hpaste(groupNames));
  t <- table(groupNames);
  verbose && cat(verbose, "Number of unique group names: ", length(t));
  verbose && print(verbose, t);
  # Sanity check
  stopifnot(length(groupNames) == nbrOfSeqs);

  keep <- (groupNames == groupName);
  bpmapList <- bpmapList[keep];
  nbrOfSeqs <- length(bpmapList);
  verbose && cat(verbose, "Number of sequences dropped: ", nbrOfSeqs0-nbrOfSeqs);
  verbose && cat(verbose, "Number of sequences after: ", nbrOfSeqs);
  verbose && cat(verbose, "Number of probes after: ", sum(getNumberOfProbes(bpmapList)));

  groupNames <- getGroupNames(bpmapList);
  verbose && cat(verbose, "Group names: ", hpaste(groupNames));
  t <- table(groupNames);
  verbose && cat(verbose, "Number of unique group names: ", length(t));
  verbose && print(verbose, t);

  verbose && exit(verbose);


  verbose && enter(verbose, "Excluding sequences that appears to be non-genomic control sequences (=all zero start positions)");
  verbose && cat(verbose, "Flavor: ", flavor);

  nbrOfSeqs0 <- nbrOfSeqs;
  verbose && cat(verbose, "Number of sequences before: ", nbrOfSeqs);
  verbose && cat(verbose, "Number of probes before: ", sum(getNumberOfProbes(bpmapList)));

  spList <- getStartPositions(bpmapList);
  verbose && cat(verbose, "Start positions:");
  verbose && str(verbose, spList);

  if (flavor == "v2") {
    keep <- sapply(spList, FUN=function(pos) !all(pos == 0L));
  } else if (flavor == "v1") {
    keep <- sapply(spList, FUN=function(pos) all(pos > 0L));
  }
  bpmapList <- bpmapList[keep];
  nbrOfSeqs <- length(bpmapList);
  verbose && cat(verbose, "Number of sequences dropped: ", nbrOfSeqs0-nbrOfSeqs);
  verbose && cat(verbose, "Number of sequences after: ", nbrOfSeqs);
  verbose && cat(verbose, "Number of probes after: ", sum(getNumberOfProbes(bpmapList)));

  groupNames <- getGroupNames(bpmapList);
  verbose && cat(verbose, "Group names: ", hpaste(groupNames));
  t <- table(groupNames);
  verbose && cat(verbose, "Number of unique group names: ", length(t));
  verbose && print(verbose, t);

  verbose && exit(verbose);

  verbose && exit(verbose);




  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Coercing to CDF friendly structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Coercing data needed for CDF structure");
  verbose && cat(verbose, "Number of sequences: ", nbrOfSeqs);
  bpmapdfList <- lapply(bpmapList, FUN=bpmapUnit2df);
  # Sanity check
  stopifnot(length(bpmapdfList) == nbrOfSeqs);
  verbose && exit(verbose);
  # Not needed anymore
  bpmapList <- NULL; # Not needed anymore



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Creating CDF tree structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating CDF tree structure for ", length(bpmapdfList), " BPMAP sequences");

  verbose && cat(verbose, "Maximum probe distance: ", maxProbeDistance);
  verbose && cat(verbose, "Minimum number of probes: ", minNbrOfProbes);
  if (!is.null(stringRemove)) {
    verbose && cat(verbose, "Regular expression used to infer unitname prefix: ", stringRemove);
  }

  # Allocate
  cdfList <- list();
  unitNames <- character(0L);

  startPositionList <- list();
  unitPrefixes <- character(0L);

  # All CDF unit will have a single group
  unitGroups <- vector("list", length=1L);

  uu <- 1L;
  for (ii in seq_along(bpmapdfList)) {
    name <- names(bpmapdfList)[ii];
    verbose && enter(verbose, sprintf("Sequence #%d ('%s') of %d", ii, name, length(bpmapdfList)));

    # Access ones
    bpmapdf <- bpmapdfList[[ii]];

    chrII <- bpmapdf$seqname[1];
    if (!is.null(stringRemove)) {
      chrII <- gsub(stringRemove, "", chrII);
    }
    verbose && cat(verbose, "Unit name prefix: ", chrII);

    nbrOfProbesII <- nrow(bpmapdf);
    verbose && cat(verbose, "Number of probes in sequence: ", nbrOfProbesII);
    nbrOfProbesII0 <- nbrOfProbesII;

    sp <- bpmapdf$startpos;
    # Sanity check
    stopifnot(all(is.finite(sp)));

    # Splitting when distances between neighboring probes are too large
    d <- diff(sp);
    keep <- (d > maxProbeDistance);
    idxsII <- which(keep);
    nbrOfSplitsII <- length(idxsII);
    if (nbrOfSplitsII > 0L) {
      verbose && printf(verbose, "Splitting into %d subsequence because there were %d too large (>%d) gaps between neighboring probes.\n", nbrOfSplitsII+1L, nbrOfSplitsII, maxProbeDistance);
    }

    # (start,end):s of subsequences
    starts <- c(1L, idxsII+1L);
    ends <- c(idxsII, nbrOfProbesII);
    counts <- (ends-starts)+1L;
##    verbose && print(verbose, data.frame(start=starts, end=ends, nbrOfProbes=counts));
    # Sanity check
    stopifnot(all(is.finite(starts)));
    stopifnot(all(is.finite(ends)));
    stopifnot(all(counts > 0L));

    # Dropping probesets with too few probes
    if (flavor == "v2") {
      keep <- (counts >= minNbrOfProbes);
    } else if (flavor == "v1") {
      keep <- (counts >= minNbrOfProbes + 2L);
    }
    rowsII <- which(keep);
    nbrOfUnitsII <- length(rowsII);
    nbrOfDroppedSeqs <- length(starts) - nbrOfUnitsII;
    if (nbrOfDroppedSeqs > 0L) {
      verbose && printf(verbose, "Dropped %d subsequences with too few (<%d) probes.\n", nbrOfDroppedSeqs, minNbrOfProbes);
    }

    verbose && cat(verbose, "Number of subsequences (=CDF units) extracted: ", nbrOfUnitsII);
    nbrOfProbesII <- sum(counts[rowsII]);
    verbose && printf(verbose, "Number of probes extracted: %d (%.2f%%) of %d\n", nbrOfProbesII, 100*nbrOfProbesII/nbrOfProbesII0, nbrOfProbesII0);

    # Access once
    pmx <- bpmapdf$pmx;
    pmy <- bpmapdf$pmy;
    mmx <- bpmapdf$mmx;
    mmy <- bpmapdf$mmy;
    isPmOnly <- all(mmx == 0L);

    # Only for verbose output
    unitNamesII <- character(length=nbrOfUnitsII);

    for (jj in seq_len(nbrOfUnitsII)) {
      rowsJJ <- rowsII[jj];
      startJJ <- starts[rowsJJ];
      endJJ <- ends[rowsJJ];

      idxsJJ <- startJJ:endJJ;
      nbrOfProbesJJ <- length(idxsJJ);
      # Sanity check
      stopifnot(nbrOfProbesJJ >= minNbrOfProbes);

      spJJ <- sp[idxsJJ];

      atom <- 0L:(nbrOfProbesJJ-1L);
      indexpos <- 0L:(nbrOfProbesJJ-1L);

      if (isPmOnly) {
        # PM only
        unitGroups[[1L]] <- list(x=pmx[idxsJJ], y=pmy[idxsJJ], pbase=rep("A", times=nbrOfProbesJJ), tbase=rep("T", times=nbrOfProbesJJ), atom=atom, indexpos=indexpos, groupdirection="sense", natoms=nbrOfProbesJJ, ncellsperatom=1L);
      } else {
        # PM+MM
        unitGroups[[1L]] <- list(x=c(pmx[idxsJJ],mmx[idxsJJ]), y=c(pmy[idxsJJ],mmy[idxsJJ]), pbase=rep("A", times=2*nbrOfProbesJJ), tbase=rep(c("T","A"), each=nbrOfProbesJJ), atom=rep(atom, times=2), indexpos=rep(indexpos, times=2), groupdirection="sense", natoms=nbrOfProbesJJ, ncellsperatom=2L);
      }

      groupNames <- sprintf("%sFROM%sTO%s", chrII, spJJ[1L], spJJ[nbrOfProbesJJ]);
      names(unitGroups) <- groupNames;

      nbrOfAtomsJJ <- sum(sapply(unitGroups, FUN=function(u) u$natoms));
      nbrOfCellsJJ <- sum(sapply(unitGroups, FUN=function(u) u$natoms*u$ncellsperatom));

      if (flavor == "v2") {
        unitIdx <- uu;
      } else if (flavor == "v1") {
        unitIdx <- ii;
      }

      cdfList[[uu]] <- list(unittype=1L, unitdirection=1L, groups=unitGroups, natoms=nbrOfAtomsJJ, ncells=nbrOfCellsJJ, ncellsperatom=nbrOfCellsJJ/nbrOfAtomsJJ, unitnumber=unitIdx);
      startPositionList[[uu]] <- spJJ;

      unitName <- groupNames;
      unitNames[uu] <- unitName;
      unitPrefixes[uu] <- chrII;

      unitNamesII[jj] <- unitName;

      # Next CDF unit
      uu <- uu + 1L;
    } # for (jj ...)
    verbose && printf(verbose, "CDF units included: [%d] %s\n", length(unitNamesII), hpaste(unitNamesII));

    verbose && exit(verbose);
  } # for (ii ...)
  verbose && exit(verbose);
  nbrOfUnits <- length(cdfList);

  # Sanity check
  stopifnot(nbrOfUnits == length(unitNames));
  stopifnot(nbrOfUnits == length(startPositionList));

  names(cdfList) <- unitNames;
  names(startPositionList) <- unitPrefixes;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Writing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Writing PPS file");
  verbose && cat(verbose, "Output pathname: ", ppsPathname);
  saveObject(startPositionList, file=ppsPathname);
  # Not needed anymore
  startPositionList <- NULL; # Not needed anymore
  verbose && exit(verbose);

  verbose && enter(verbose, "Writing CDF file");
  verbose && cat(verbose, "Output pathname: ", cdfPathname);

  cdfHeader <- list(probesets=nbrOfUnits, qcprobesets=0L, reference="", chiptype=chipType, filename=cdfPathname, nqcunits=0L, nunits=nbrOfUnits, rows=rows, cols=cols, refseq="", nrows=rows, ncols=cols);
  verbose && cat(verbose, "CDF header:");
  verbose && str(verbose, cdfHeader);

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(cdfPathname, verbose=verbose);

  .writeCdf(pathnameT, cdfheader=cdfHeader, cdf=cdfList, cdfqc=NULL, overwrite=TRUE, verbose=verbose);

  # Rename temporary file
  popTemporaryFile(pathnameT, verbose=verbose);

  ## Create checksum file
  dfZ <- getChecksumFile(cdfPathname)

  verbose && exit(verbose);
  verbose && exit(verbose);

  invisible(pathname);
}) # bpmapCluster2Cdf()


############################################################################
# HISTORY:
# 2012-06-21 [HB]
# o GENERALIZATION: Now the elements of the written start position list
#   (PPS) have names corresponding to the unit prefix (e.g. "chrX").
# o Added argument 'path' to bpmapCluster2Cdf(), which now defaults to
#   annotationData/chipTypes/<chipType>/.
# o Added more internal sanity checks.
# o CLARIFICATION: Restructured the bpmapCluster2Cdf() method such that
#   is more clear how BPMAP sequences are filtered out, i.e. keeping
#   sequencing with a matching group name and excluding those that
#   appears to be non-genomic control sequences.
# o ROBUSTNESS: Arguments 'rows' and 'cols' for bpmapCluster2Cdf() are
#   mandatory (again).  The reason for this is that the BPMAP file is only
#   useful to infer a lower bound for them, but not their exact values.
# o BUG FIX: bpmapCluster2Cdf(..., minNbrOfProbes=n) filtered out units
#   with less than (n+2L) probes, not n probes.
# o BUG FIX: Previously non-genomic control sequences, which were filtered
#   out, were identified as having at least one probe start position to
#   be zero.  Now they are indentified by all start positions being zero.
#   This change was discussed today with MR via email.
# o BUG FIX: The generated CDF structure had 'unitnumber':s set to be
#   equal to BPMAP sequence index rather than the CDF unit number.  This
#   probably didn't matter for the written CDF, because writeCdf() sets
#   the unit indices itself.
# 2011-08-31 [HB]
# o All arguments after '...' except 'verbose' may be dropped in a future
#   release, especially 'rows' and 'cols'.
# o CLARIFICATION: Renamed argument 'nProbes' to 'minNbrOfProbes'.
# 2011-08-30 [HB]
# o Now bpmapCluster2Cdf() returns the pathname to the created CDF.
# 2011-08-29 [HB]
# o Replaced argument 'cdfName' with 'chipType' and 'tags'.  This
#   simplifies adding custom tags to the CDF.
#   It will also allow us (later) to write the CDF to the correct
#   annotationData/ directory.
# o ROBUSTNESS: Now the writing of the CDF file is atomic by first writing
#   to a temporary file which is then renamed.
# o GENERALIZATION/ROBUSTNESS: Now the method precount the number of
#   CDF units so it can allocate a CDF tree structure of the correct size.
# o GENERALIZATION: Now the default value of argument 'groupName' is
#   more general.  Before it was hardcoded to "Hs".
# o GENERALIZATION: Now the default value of argument 'stringRemove' is
#   more general.  Before it was hardcoded to a particular BPMAP file.
# o DOCUMENTATION: Added some basic Rdoc documentation.
# o Added more verbose progress output.
# o CLEANUP: Allocated variable 'll' was never used.
# 2011-07-11 [HB]
# o CLEANUP: Code and argument cleanup.
# o ROBUSTNESS: Using '&&' instead of '&' if-statement.
# o ROBUSTNESS: Added more argument validation.
# o ROBUSTNESS: Added protection against overwriting existing CDFs.
# o ROBUSTNESS: Added tests for valid CDF filenames.
# 2008-11-28 [HB]
# o CLEANUP: Tidying up code.
# o Now using R.utils::Verbose statements.
# 2008-11-xx [MR]
# o Created.
############################################################################
