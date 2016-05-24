###########################################################################/**
# @set "class=ChipEffectFile"
# @RdocMethod getAM
#
# @title "Gets the log-intensities and log-ratios of chip effects for two arrays"
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
#  Returns a Nx2 matrix where N is the number of units returned.
#  The names of the columns are A (log-intensities) and M (log-ratios).
#  The names of the rows are the unit indices (as indexed by the CDF).
#  The rows are ordered according to \code{units} arguments.
# }
#
# @author "HB"
#
# \seealso{
#   @seemethod "getXAM".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getAM", "ChipEffectFile", function(this, other, units=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  cdf <- getCdf(this);
  # Argument 'other':
  if (is.null(other)) {
    # Do not calculate ratios relative to a reference
  } else {
    other <- Arguments$getInstanceOf(other, "ChipEffectFile");
  }

  # Argument 'units':
  ugcMap <- NULL;
  if (is.null(units)) {
  } else if (isUnitGroupCellMap(units)) {
    ugcMap <- units;
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits(cdf));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Getting (A,M)-transformed chip effects");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get cell map
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (is.null(ugcMap)) {
    verbose && enter(verbose, "Get unit-to-cell map");
    ugcMap <- getUnitGroupCellMap(this, units=units, verbose=less(verbose));
    verbose && exit(verbose);
  }
  # Not needed anymore
  units <- NULL;
  nunits <- nrow(ugcMap);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the sample
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving sample thetas");
  theta <- getDataFlat(this, units=ugcMap, fields="theta",
                                          verbose=less(verbose))[,"theta"];
  nTheta <- length(theta);
  if (!identical(nTheta, nunits)) {
    verbose && str(verbose, theta);
    verbose && print(verbose, nunits);
    throw("Internal error: The number of chip-effect values is not equal to the number of units requested: ", nTheta, " != ", nunits);
  }
  verbose && exit(verbose);

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
    thetaR <- getDataFlat(other, units=ugcMap, fields="theta",
                                            verbose=less(verbose))[,"theta"];
    stopifnot(identical(length(thetaR), nTheta));
    verbose && exit(verbose);
  } # if (!is.null(other))

  # Calculate raw copy numbers relative to reference?
  # AD HOC /HB 2009-11-22 (get[X]AM() should be dropped in the future)
  if (is.null(other)) {
    M <- theta;
    A <- theta;
  } else {
    M <- theta / thetaR;
    A <- theta * thetaR;
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Log2 scale
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  M <- log(M, base=2);
  A <- log(A, base=2)/2;
  stopifnot(identical(length(M), nTheta));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the return matrix
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  am <- matrix(c(A,M), ncol=2)
  colnames(am) <- c("A", "M");
  rownames(am) <- ugcMap[,"unit"];

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
setMethodS3("getXAM", "ChipEffectFile", function(this, other, chromosome, units=NULL, ..., verbose=FALSE) {
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
  chromosome <- Arguments$getCharacter(chromosome);  # integer? /HB 2009-05-18

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
  x <- getPositions(gi, units=units, verbose=less(verbose));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove units for which we have no position information
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
# 2007-03-15
# o Replaced "SNPs" with "units" in error messages.
# 2007-01-16
# o BUG FIX: getAM() and getXAM() required that argument 'other' was an
#   CnChipEffectFile; ChipEffectFile is enough.
# 2006-12-02
# o BUG FIX: getAM() did not handle units=NULL.
# 2006-11-28
# o Now getAM() is making use of the new getCellMap() and getDataFlat()
#   functions to speed up the reading.
# 2006-10-31
# o Added Rdoc comments.
# 2006-09-26
# o Created.
############################################################################
