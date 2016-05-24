###########################################################################/**
# @RdocClass HetLogAddPlm
#
# @title "The HetLogAddPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents a log-additive model similar to the one described
#  in Irizarry et al (2003), except that the errors may have different
#  variances for different probes.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaPlm".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#  @see "RmaPlm".
# }
#*/###########################################################################
setConstructorS3("HetLogAddPlm", function(...) {
  extend(RmaPlm(...), "HetLogAddPlm")
})


setMethodS3("getAsteriskTags", "HetLogAddPlm", function(this, collapse=NULL, ...) {
  # Returns 'RMA[,<flavor>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "HLA";

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getFitUnitGroupFunction
#
# @title "Gets the low-level function that fits the PLM"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @function.
# }
#
# \author{
#   Henrik Bengtsson and Ken Simpson (WEHI) utilizing Ben Bolstad's
#   \pkg{preprocessCore} package.
# }
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitUnitGroupFunction", "HetLogAddPlm", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting the PLM fit function");

  # Early error, if package is missing
  requireNamespace("preprocessCore") || throw("Package not loaded: preprocessCore")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the flavor of fitting algorithm for the RMA PLM
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Shift signals?
  shift <- this$shift;
  if (is.null(shift)) {
    shift <- 0;
  }
  verbose && cat(verbose, "Amount of shift: ", shift);

  fitPlm <- function(y, ...) {
    fit <- fitWHLAPLM.matrix(y+shift, ...);
    fit$sdTheta <- fit$seTheta;
    fit$sdPhi <- fit$sePhi;
    fit$thetaOutliers <- rep(NA_real_, times=ncol(y));
    fit$phiOutliers <- rep(NA_real_, times=nrow(y));
    if (length(fit$sdTheta) != ncol(y)) {
      print(list(y=y, fit=fit));
      throw("Internal error");
    }
    if (length(fit$sdPhi) != nrow(y)) {
      print(list(y=y, fit=fit));
      throw("Internal error");
    }
    fit;
  }


  verbose && str(verbose, fitPlm);

  # Test that it works and is available.
  verbose && enter(verbose, "Validating the fit function on some dummy data");
  ok <- FALSE;
  tryCatch({
    fitPlm(matrix(1:6+0.1, ncol=3));
    ok <- TRUE;
  }, error = function(ex) {
    print(ex);
  })
  if (!ok) {
    throw("The fit function failed");
  }
  verbose && exit(verbose);


  verbose && exit(verbose);

  fitPlm;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-08-21
# o Now getFitUnitGroupFunction() for HetLogAddPlm allocates using
#   as.double(NA) instead of NA.
# 2007-10-05
# o Created from RmaPlm.R.
############################################################################
