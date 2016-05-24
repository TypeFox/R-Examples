###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod getAM
#
# @title "Gets the log-intensities and log-ratios of chip effects of the set relative to a reference chip effect file"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{other}{The second @see "ChipEffectFile" object used as the
#     reference.}
#   \item{units}{(The subset of units to be matched.
#     If @NULL, all units are considered.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an Jx2xK @array where J is the number of units, and K is
#  the number of arrays (arrays are always the last dimension).
#  The names of the columns are A (log-intensities) and M (log-ratios).
#  The names of the rows are the unit indices (as indexed by the CDF).
#  The rows are ordered according to \code{units} arguments.
# }
#
# @author "HB"
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAM", "ChipEffectSet", function(this, other, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  if (is.null(other)) {
    # Do not calculate ratios relative to a reference
  } else {
    other <- Arguments$getInstanceOf(other, "ChipEffectFile");
  }

  # Argument 'units':
  cdf <- getCdf(this);
  if (is.null(units)) {
    nunits <- nbrOfUnits(cdf);
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
    nunits <- length(units);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting (A,M)-transformed chip effects for a set of arrays");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get cell map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Get unit-to-cell map");
  cf <- getOneFile(this, mustExist=TRUE);
  map <- getUnitGroupCellMap(cf, units=units, verbose=less(verbose));
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate return array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  arrayNames <- getNames(this);
  dimnames <- list(map[,"unit"], arrayNames, c("A", "M"));
  dim <- sapply(dimnames, FUN=length);
  am <- array(NA, dim=dim, dimnames=dimnames);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the other
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(other)) {
    verbose && enter(verbose, "Retrieving other thetas");

    # Workaround for now (just in case). /HB 2006-09-26 TODO
    if (inherits(other, "SnpChipEffectFile")) {
      other$mergeStrands <- this$mergeStrands;
      if (inherits(other, "CnChipEffectFile")) {
        other$combineAlleles <- this$combineAlleles;
      }
    }

    # Get the other theta estimates
    thetaR <- getDataFlat(other, units=map, fields="theta", verbose=less(verbose))[,"theta"];
    nTheta <- length(thetaR);
    stopifnot(identical(nTheta, nrow(map)));
    verbose && exit(verbose);
  } # if (!is.null(other))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving sample thetas");
  for (aa in seq_along(this)) {
    cf <- this[[aa]];
    theta <- getDataFlat(cf, units=map, fields="theta", verbose=less(verbose))[,"theta"];
    if (!identical(length(theta), nTheta)) {
      verbose && str(verbose, theta);
      verbose && print(verbose, nunits);
      throw("Internal error: The number of chip-effect values is not equal to the number of units requested: ", length(theta), " != ", nTheta);
    }

    # Calculate raw copy numbers relative to reference?
    # AD HOC /HB 2009-11-22 (get[X]AM() should be dropped in the future)
    if (is.null(other)) {
      M <- theta;
      A <- theta;
    } else {
      M <- theta / thetaR;
      A <- theta * thetaR;
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Log2 scale
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    M <- log(M, base=2);
    A <- log(A, base=2)/2;
    stopifnot(identical(length(M), nTheta));

    am[,aa,"A"] <- A;
    am[,aa,"M"] <- M;

    # Not needed anymore
    theta <- A <- M <- NULL;
  } # for (aa in ...)
  verbose && exit(verbose);

  verbose && exit(verbose);

  am;
}) # getAM()



###########################################################################/**
# @RdocMethod getXAM
#
# @title "Gets the physical position, log-intensities and log-ratios of chip effects for two arrays"
#
# \description{
#  @get "title" of units on a certain chromosome.
# }
#
# @synopsis
#
# \arguments{
#   \item{other}{The second @see "ChipEffectFile" object used as the
#     reference.}
#   \item{chromosome}{(The chromosome for which results should be returned.}
#   \item{units}{(The subset of units to be matched.
#     If @NULL, all units are considered.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a Nx3 matrix where N is the number of units returned.
#  The names of the columns are X (physical position in a given chromosome),
#  A (log-intensities) and M (log-ratios).
#  The names of the rows are the unit indices (as indexed by the CDF).
#  \emph{Note: The rows are ordered according to chromosomal position.}
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "getAM".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getXAM", "ChipEffectSet", function(this, other, chromosome, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'other':
  if (is.null(other)) {
    # Do not calculate ratios relative to a reference
  } else {
    other <- Arguments$getInstanceOf(other, "ChipEffectFile");
  }

  # Argument 'chromosome':
  chromosome <- Arguments$getChromosome(chromosome);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting (X,A,M)-transformed chip effects");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieve genome information, i.e. chromosome positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving genome information");
  cdf <- getCdf(this);
  gi <- getGenomeInformation(cdf);

  # Extract the units from the given chromosome.  Requested units not on
  # chromosome are ignored.
  if (!is.null(units))
    verbose && str(verbose, units);
  units <- getUnitIndices(gi, chromosome=chromosome, units=units, verbose=less(verbose));
  nunits <- length(units);
  if (nunits == 0)
    throw("No units found on requested chromosome: ", chromosome);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the relative copy-number estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  am <- getAM(this, other=other, units=units, verbose=less(verbose));

  # Get the unit indices for all unit groups
  units <- as.integer(rownames(am));

  # Get the positions of all unit groups
  x <- getPositions(gi, units=units);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove SNPs for which we have no position information
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keep <- which(!is.na(x));
  nexcl <- length(x) - length(keep);
  if (nexcl > 0) {
    msg <- sprintf("Could not find position information on %d unit groups: ", nexcl);
    verbose && cat(verbose, msg);
    warning(msg);
    x <- x[keep];
    units <- units[keep];
  }

  am <- am[,c("M","A"), drop=FALSE]; # Ad hoc /HB 2007-02-19
  xam <- cbind(x=x, am);

  verbose && exit(verbose);

  xam;
}, protected=TRUE) # getXAM()


############################################################################
# HISTORY:
# 2009-11-22
# o Now get[X]AM() accepts other=NULL.
# 2007-03-04
# o Added getAM().
# o Created from ChipEffectFile.xam.R.
############################################################################
