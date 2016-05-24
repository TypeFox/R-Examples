###########################################################################/**
# @RdocClass ExonRmaPlm
#
# @title "The ExonRmaPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model part of the Robust Multichip
#  Analysis (RMA) method described in Irizarry et al (2003), as implemented
#  for exon arrays.  The model may be fitted with exons merged into
#  transcripts (all probes fitted together) or on an individual exon basis
#  (probes within an exon treated as a group, but exons fitted separately).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "RmaPlm".}
#   \item{mergeGroups}{A @logical flag specifying whether to merge exons
#      into transcripts.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#    @see "RmaPlm".
# }
#
# @author "KS, HB, EP"
#
# \references{
#  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}.
#  NAR, 2003, 31, e15.\cr
# }
#*/###########################################################################
setConstructorS3("ExonRmaPlm", function(..., mergeGroups=TRUE) {
  extend(RmaPlm(...), "ExonRmaPlm",
    mergeGroups=mergeGroups
  )
})



setMethodS3("getAsteriskTags", "ExonRmaPlm", function(this, collapse=NULL, ...) {
  # Returns 'RMA[,<flavor>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);

  # Add class specific parameter tags
  if (this$mergeGroups)
    tags <- c(tags, "merged");

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


# utility function - keep here for now

cdfMergeGroups <- function(groups, ...) {
  nbrOfGroups <- length(groups);

  nbrOfFields <- length(.subset2(groups,1));
  newGroup <- vector("list", nbrOfFields);
  for (ff in seq_len(nbrOfFields)) {
    newGroup[[ff]] <- unlist(lapply(groups, FUN=.subset2, ff), use.names=FALSE);
  }
  names(newGroup) <- names(.subset2(groups,1));
  return(list(newGroup));
}

setMethodS3("getCellIndices", "ExonRmaPlm", function(this, ...) {

  cells <- NextMethod("getCellIndices");

  # Merge groups?
  if (this$mergeGroups) {
    cells <- .applyCdfGroups(cells, cdfMergeGroups);
  }

  cells;
})


setMethodS3("getChipEffectSetClass", "ExonRmaPlm", function(this, ...) {
  ExonChipEffectSet;
}, private=TRUE)


setMethodS3("getChipEffectSet", "ExonRmaPlm", function(this, ...) {
  ces <- NextMethod("getChipEffectSet");
  setMergeGroups(ces, this$mergeGroups);
  ces;
})


setMethodS3("getProbeAffinityFile", "ExonRmaPlm", function(this, ..., .class=ExonProbeAffinityFile) {
  paf <- NextMethod("getProbeAffinityFile", .class=.class);
  setMergeGroups(paf, this$mergeGroups);
  paf;
})


setMethodS3("setMergeGroups", "ExonRmaPlm", function(this, ...) {
  ces <- getChipEffectSet(this);
  setMergeGroups(ces, ...);
  paf <- getProbeAffinityFile(this);
  setMergeGroups(paf, ...);
})



setMethodS3("getParameters", "ExonRmaPlm", function(this, ...) {
  params <- NextMethod("getParameters");
  params$mergeGroups <- this$mergeGroups;
  params;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod getFitUnitGroupFunction
#
# @title "Gets the low-level function that fits the Exon PLM"
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
setMethodS3("getFitUnitGroupFunction", "ExonRmaPlm", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Thresholds for skipping/using median polish
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  skipThreshold <- getOption(aromaSettings,
                                 "models/RmaPlm/skipThreshold", c(Inf, Inf));

  medianPolishThreshold <- getOption(aromaSettings,
                         "models/RmaPlm/medianPolishThreshold", c(Inf, Inf));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rmaModelAffyPlm()
  # Author: Henrik Bengtsson, UC Berkeley.
  # Requires: affyPLM() by Ben Bolstad.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  exonRmaModel <- function(y, psiCode=0, psiK=1.345, ...){
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Assert right dimensions of 'y'.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # If input data are dimensionless, return NAs. /KS 2006-01-30
    dim <- dim(y);
    if (is.null(dim)) {
      nbrOfArrays <- length(getDataSet(this));
      return(list(theta=rep(NA, nbrOfArrays),
                  sdTheta=rep(NA, nbrOfArrays),
                  thetaOutliers=rep(NA, nbrOfArrays),
                  phi=c(),
                  sdPhi=c(),
                  phiOutliers=c()));

    }

    if (length(dim) != 2) {
      stop("Argument 'y' must have two dimensions: ",
                                                paste(dim, collapse="x"));
    }

    K <- dim[1];  # Number of probes
    I <- dim[2];  # Number of arrays

    # Skip unit group?
    if (K > skipThreshold[1] && I > skipThreshold[2]) {
      warning("Ignoring a unit group when fitting probe-level model, because it has a ridiculously large number of data points: ", paste(dim, collapse="x"), " > ", paste(skipThreshold, collapse="x"));

      return(list(theta=rep(NA, I),
                  sdTheta=rep(NA, I),
                  thetaOutliers=rep(NA, I),
                  phi=rep(NA, K),
                  sdPhi=rep(NA, K),
                  phiOutliers=rep(NA, K)
                 )
            );
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Transform data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Add shift?
    if (shift != 0)
      y <- y + shift;

    # Log-additive model
    y <- log(y, base=2);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Use median polish for large probesets?
    if (K > medianPolishThreshold[1] && I > medianPolishThreshold[2]) {
      mp <- medpolish(y, trace.iter=FALSE);
      fit <- list(
        Estimates = c(mp$overall+mp$col, mp$row),
        StdErrors = rep(0, length(c(mp$row, mp$col)))
      );
    } else {
      # Fit model using affyPLM code
      fit <- rlm(y, psiCode, psiK)
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extract parameters
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    est <- fit$Estimates;
    se <- fit$StdErrors;

    # Chip effects
    beta <- est[1:I];

    # Probe affinities.  If only one probe, must have affinity=1 since
    # sum constraint => affinities sum to zero (on log scale)
    if (K == 1) {
      alpha <- 0;
    } else {
      alpha <- est[(I+1):length(est)];
      alpha[length(alpha)] <- -sum(alpha[1:(length(alpha)-1)]);
    }

    # Estimates on the intensity scale
    theta <- 2^beta;
    phi <- 2^alpha;

    # The RMA model is fitted with constraint sum(alpha) = 0, that is,
    # such that prod(phi) = 1.

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fit function must return: theta, sdTheta, thetaOutliers,
    # phi, sdPhi, phiOutliers.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(se)) {
      # For affyPLM v1.10.0 (2006-09-26) or older.
      sdTheta <- rep(1, I);
      sdPhi <- rep(1, K);
    } else {
      # For affyPLM v1.11.6 (2006-11-01) or newer.
      sdTheta <- 2^(se[1:I]);
      sdPhi <- 2^(se[(I+1):length(se)]);
    }
    thetaOutliers <- rep(FALSE, I);
    phiOutliers <- rep(FALSE, K);

    # Return data on the intensity scale
    list(theta=theta, sdTheta=sdTheta, thetaOutliers=thetaOutliers,
         phi=phi, sdPhi=sdPhi, phiOutliers=phiOutliers);
  } # rmaModelAffyPlm()
  attr(exonRmaModel, "name") <- "exonRmaModel";



  getRlmPkg <- function(..., verbose=FALSE) {
    # First, try to see if preprocessCore > v0.99.14 is available
    pkg <- "preprocessCore";
    pkgDesc <- packageDescription(pkg);
    if (is.list(pkgDesc)) {
      ver <- pkgDesc$Version;
      verbose && cat(verbose, pkg, " version: ", ver);
      if (compareVersion(ver, "0.99.14") >= 0)
        return(pkg);
    }

    # Second, try to see if affyPLM <= v1.13.8 is available
    pkg <- "affyPLM";
    pkgDesc <- packageDescription(pkg);
    if (is.list(pkgDesc)) {
      ver <- pkgDesc$Version;
      verbose && cat(verbose, pkg, " version: ", ver);
      if (compareVersion(ver, "1.13.8") <= 0)
        return(pkg);
    }

    throw("Neither preprocessCore v0.99.14+ nor affyPLM v1.13.8- is available.");
  } # getRlmPkg()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting the PLM fit function");

  # Shift signals?
  shift <- this$shift;
  if (is.null(shift))
    shift <- 0;
  verbose && cat(verbose, "Amount of shift: ", shift);


  fcnList <- RmaPlm$getRlmFitFunctions(verbose=less(verbose));
  verbose && str(verbose, fcnList);
  # To please R CMD check
  rlm <- wrlm <- NULL; rm(list=c("rlm", "wrlm"));
  attachLocally(fcnList);
  rmaModel <- exonRmaModel;
  verbose && str(verbose, rmaModel);

  # Test that it works and is available.
  verbose && enter(verbose, "Validating the fit function on some dummy data");
  ok <- FALSE;
  tryCatch({
    rmaModel(matrix(1:6+0.1, ncol=3));
    ok <- TRUE;
  }, error = function(ex) {
    print(ex);
  })
  if (!ok) {
    throw("The fit function for requested exon RMA PLM failed");
  }
  verbose && exit(verbose);

  verbose && exit(verbose);

  rmaModel;
}, private=TRUE)




##############################################################################
# HISTORY:
# 2008-06-03 [HB]
# o Added getParameterSet() to ExonRmaPlm so that 'mergeGroups' is returned.
# 2007-12-02
# o Added mechanism to fit "large" unit groups with median polish (for speed).
# o Added mechanism to avoid fitting unit groups with ridiculously many cells.
# 2007-12-06
# o Added getAsteriskTag() for ExonRmaPlm.
# 2007-11-07
# o Now the fit function supports 'shift' too. Before it was ignored.
# 2007-09-20
# o Updated getFitFunction() for ExonRmaPlm to deal with compatibility issues
#   related different versions of affyPLM/preprocessCore.
# o Added verbose output to getFitFunction().
# o Renamed the variables such that index I is for samples and K is for
#   probes, as in the paper.
# 2007-07-13 /HB
# o Removed findUnitsTodo() of ExonRmaPlm; it just replicated the superclass.
# 2007-04-24 /HB+EP
# o getProbeAffinityFile() of ExonRmaPlm did not return the correct subclass.
# 2007-04-13 /HB+EP
# o BUG FIX: getChipEffectSet() and getProbeAffinityFile() did not set the
#   'mergeStrands' parameter.  Thanks Elizabeth Purdom for the fix.
# 2006-??-??
# o Created by Ken Simpson, WEHI.
##############################################################################
