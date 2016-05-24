########################################################################/**
# @RdocDefault gtypeCelToPQ
#
# @title "Function to immitate Affymetrix' gtype\_cel\_to\_pq software"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{filename}{The name of a CEL file.}
#   \item{units}{Indices of CDF units to be returned.}
#   \item{...}{Arguments passed to @see "affxparser::readCelUnits".}
#   \item{cdf}{A CDF @list structure, the pathname of the CDF file, or
#     @NULL.  If @NULL, the CDF file corresponding to the chip type of
#     the CEL file is searched for using @see "affxparser::findCdf".}
#   \item{nbrOfQuartets}{The number of probe quartets in the returned
#     @matrix.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an NxK @matrix where N is the number of probesets (SNPs) and
#  K=4*Q where Q is the number of probe quartets (PMA,MMA,PMB,MMB).
#  The rownames corresponds to the probeset names.
# }
#
# @examples "../incl/gtypeCelToPQ.Rex"
#
# \seealso{
#  @see "gtypeCelToPQ".
#  @see "affxparser::applyCdfGroups".
# }
#
# @author
#
# \references{
#  [1] Affymetrix, \emph{Genotyping Probe Set Structure},
#      Developers' Network, White paper, 2005-2015.
# }
#
# @keyword programming
#**/#######################################################################
setMethodS3("gtypeCelToPQ", "default", function(filename, units=NULL, ..., cdf=NULL, nbrOfQuartets=NULL, verbose=FALSE) {
  requireNamespace("affxparser") || throw("Package not loaded: affxparser");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename':
  filename <- Arguments$getReadablePathname(filename, mustExist=TRUE);

  # Argument 'nbrOfQuartets':
  if (!is.null(nbrOfQuartets)) {
    nbrOfQuartets <- Arguments$getInteger(nbrOfQuartets, range=c(1,1024));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get CDF structure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.list(cdf)) {
  } else if (is.character(cdf)) {
    cdfFile <- cdf;
  } else {
    verbose && enter(verbose, "Obtaining chip type from CEL header");
    celHeader <- affxparser::readCelHeader(filename);
    chiptype <- celHeader$chiptype;
    verbose && exit(verbose);
    verbose && enter(verbose, "Searching for CDF file");
    cdfFile <- affxparser::findCdf(chiptype);
    verbose && exit(verbose);
  }

  if (!is.list(cdf)) {
    verbose && enter(verbose, "Reading CDF structure");
    cdf <- affxparser::readCdfUnits(cdfFile, units=units, stratifyBy="pmmm");
    verbose && exit(verbose);

    verbose && enter(verbose, "Rearranging CDF structure");
    cdfGtypeCelToPQ <- affxparser::cdfGtypeCelToPQ;
    cdf <- affxparser::applyCdfGroups(cdf, cdfGtypeCelToPQ);
    verbose && exit(verbose);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read CEL intensities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading CEL file");
  cel <- affxparser::readCelUnits(filename, ..., cdf=cdf);
  verbose && exit(verbose);

  cdf <- NULL; # Not needed anymore



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create (PMA,MMA,PMB,MMB) matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Creating matrix");
  # Assume equal number of probes in all probesets, which is also what
  # Affymetrix' gtype_cel_to_pq does.
  if (is.null(nbrOfQuartets)) {
    nbrOfProbes <- length(cel[[1]][[1]][[1]]);
    # Assume a multiple of four probes
    if (nbrOfProbes %% 4 != 0) {
      throw("The first unit does not have a multiple of four (quartet) probes: ", nbrOfProbes);
    }
    nbrOfQuartets <- nbrOfProbes %/% 4;
  }

  nbrOfUnits <- length(cel);
  rownames <- names(cel);
  quartets <- rep(c("pmA", "mmA", "pmB", "mmB"), times=nbrOfQuartets);
  colnames <- paste(quartets, rep(sprintf("%02d", 1:nbrOfQuartets), each=4), sep="");
  dimnames <- list(rownames, colnames);
  res <- matrix(NA, nrow=nbrOfUnits, ncol=4*nbrOfQuartets, dimnames=dimnames);
  cc <- 1:ncol(res);
  for (uu in 1:nbrOfUnits) {
    res[uu,cc] <- cel[[uu]][[1]][[1]][cc];
  }

  verbose && exit(verbose);

  res;
})

############################################################################
# HISTORY:
# 2006-03-12
# o Added argument 'nbrOfQuartets'.
# 2006-03-11
# o Tried it.  It can be used to generate the same tab-delimited files as
#   Affymetrix' gtype_cel_to_pq application.
# o Created, partly to verify the correctness of affxparser.
############################################################################
