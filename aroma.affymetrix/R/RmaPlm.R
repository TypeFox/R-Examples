###########################################################################/**
# @RdocClass RmaPlm
#
# @title "The RmaPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model part of the Robust Multichip
#  Analysis (RMA) method described in Irizarry et al (2003).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
#   \item{flavor}{A @character string specifying what model fitting algorithm
#     to be used.  This makes it possible to get identical estimates as other
#     packages.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   For a single unit group, the log-additive model of RMA is:
#
#    \deqn{log_2(y_{ik}) = \beta_i + \alpha_k + \varepsilon_{ik}}
#
#   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I},
#   and \eqn{\alpha_k} are the probe affinities for probes \eqn{k=1,...,K}.
#   The \eqn{\varepsilon_{ik}} are zero-mean noise with equal variance.
#   The model is constrained such that \eqn{\sum_k{\alpha_k} = 0}.
#
#   Note that all PLM classes must return parameters on the intensity scale.
#   For this class that means that \eqn{\theta_i = 2^\beta_i} and
#   \eqn{\phi_k = 2^\alpha_k} are returned.
# }
#
# \section{Different flavors of model fitting}{
#   There are a few differ algorithms available for fitting the same
#   probe-level model.  The default and recommended method
#   (\code{flavor="affyPLM"}) uses the implementation in the
#   \pkg{preprocessCore} package which fits the model parameters robustly
#   using an M-estimator (the method used to be in \pkg{affyPLM}).
#
#   Alternatively, other model-fitting algorithms are available.
#   The algorithm (\code{flavor="oligo"}) used by the \pkg{oligo} package,
#   which originates from the \pkg{affy} packages, fits the model using
#   median polish, which is a non-robust estimator.  Note that this algorithm
#   does not constraint the probe-effect parameters to multiply to one on
#   the intensity scale.  Since the internal function does not return these
#   estimates, we can neither rescale them.
# }
#
# @author "HB, KS"
#
# \references{
#  Irizarry et al. \emph{Summaries of Affymetrix GeneChip probe level data}.
#  NAR, 2003, 31, e15.\cr
# }
#*/###########################################################################
setConstructorS3("RmaPlm", function(..., flavor=c("affyPLM", "oligo")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  this <- extend(ProbeLevelModel(...), "RmaPlm",
    .flavor = flavor,
    treatNAsAs = "weights"
  );
  validate(this);
  this;
})


setMethodS3("getAsteriskTags", "RmaPlm", function(this, collapse=NULL, ...) {
  # Returns 'PLM[,<shift>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "RMA";

  # Add class specific parameter tags
  if (this$.flavor != "affyPLM")
    tags <- c(tags, this$.flavor);

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getParameters", "RmaPlm", function(this, ...) {
  params <- NextMethod("getParameters");
  params$flavor <- this$.flavor;
  params$treatNAsAs <- this$treatNAsAs;
  params;
}, protected=TRUE)



setMethodS3("getProbeAffinityFile", "RmaPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the probe affinities (and create files etc)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paf <- NextMethod("getProbeAffinityFile");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the encode and decode functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setEncodeFunction(paf, function(groupData, ...) {
    phi <- .subset2(groupData, "phi");
    stdvs <- .subset2(groupData, "sdPhi");
    outliers <- .subset2(groupData, "phiOutliers");

    # Encode outliers as the sign of 'pixels'; -1 = TRUE, +1 = FALSE
    pixels <- sign(0.5 - as.integer(outliers));

    list(intensities=phi, stdvs=stdvs, pixels=pixels);
  })

  setEncodeFunction(paf, function(groupData, ...) {
    list(
      intensities = .subset2(groupData, "phi"),
      stdvs = .subset2(groupData, "sdPhi"),
      # Encode outliers as the sign of 'pixels'; -1 = TRUE, +1 = FALSE
      pixels = ifelse(.subset2(groupData, "phiOutliers"), -1, +1)
    );
  })

##  setEncodeFunction(paf, function(groupData, ...) {
##    groupData[[3]] <- ifelse(.subset2(groupData, "phiOutliers"), -1, +1);
##    names(groupData) <- c("phi", "sdPhi", "pixels");
##    groupData;
##  })
##
##  setEncodeFunction(paf, function(groupData, ...) {
##    groupData[[3]] <- -1*.subset2(groupData, 3);
##    names(groupData) <- c("phi", "sdPhi", "pixels");
##    groupData;
##  })

  setDecodeFunction(paf,  function(groupData, ...) {
    intensities <- .subset2(groupData, "intensities");
    stdvs <- .subset2(groupData, "stdvs");
    pixels <- .subset2(groupData, "pixels");

    # Outliers are encoded by the sign of 'pixels'.
    outliers <- as.logical(1-sign(pixels));

    list(
      phi=intensities,
      sdPhi=stdvs,
      phiOutliers=outliers
    );
  })

  paf;
}, private=TRUE)



setMethodS3("getRlmFitFunctions", "RmaPlm", function(static, withPriors=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'withPriors':
  withPriors <- Arguments$getLogical(withPriors);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  fcnList <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # First, try the preprocessCore package
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pkg <- "preprocessCore";
  pkgDesc <- packageDescription(pkg);
  if (is.list(pkgDesc)) {
    ver <- pkgDesc$Version;
    verbose && cat(verbose, pkg, " version: ", ver);

    # Internal sanity check
    if (compareVersion(ver, "1.8.0") < 0) {
      throw("Requires preprocessCore >= 1.8.0");
    }

    if (withPriors) {
      # To please R CMD check
      PACKAGE <- "preprocessCore";
      fcnList <- list(
        wrlm = function(y, phi, psiCode, psiK, w, scale=NULL) {
          # From preprocessCore::rcModelWPLM()
          .Call("R_rlm_rma_given_probe_effects",
                y, phi, psiCode, psiK, w, scale, PACKAGE=PACKAGE);
        },
        rlm = function(y, phi, psiCode, psiK, w, scale=NULL) {
          # From preprocessCore::rcModelPLM()
          .Call("R_rlm_rma_given_probe_effects",
                y, phi, psiCode, psiK, scale, PACKAGE=PACKAGE);
        }
      );
    } else {
      # To please R CMD check
      PACKAGE <- "preprocessCore";
      fcnList <- list(
        wrlm = function(y, psiCode, psiK, w, scale=NULL) {
          # From preprocessCore::rcModelPLM()
          .Call("R_wrlm_rma_default_model",
                y, psiCode, psiK, w, scale, PACKAGE=PACKAGE);
        },
        rlm = function(y, psiCode, psiK, w, scale=NULL) {
          # From preprocessCore::rcModelPLM()
          .Call("R_rlm_rma_default_model",
                y, psiCode, psiK, scale, PACKAGE=PACKAGE);
        }
      );
    }
  }

  if (!is.null(fcnList)) {
    .require <- require
    .require(pkg, character.only=TRUE) || throw("Package not loaded: ", pkg)
    return(fcnList);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Failure
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  throw("Failed to retrieve RLM fit functions. Please install the 'preprocessCore' package.");
}, static=TRUE, protected=TRUE) # getRlmFitFunctions()



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
setMethodS3("getFitUnitGroupFunction", "RmaPlm", function(this, ..., verbose=FALSE) {
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

  flavor <- this$.flavor;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rmaModelAffyPlm()
  # Author: Henrik Bengtsson, UC Berkeley.
  # Requires: affyPLM() by Ben Bolstad.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rmaModelAffyPlm <- function(y, priors=NULL, ..., psiCode=0, psiK=1.345){
    # Assert right dimensions of 'y'.

    # If input data are dimensionless, return NAs. /KS 2006-01-30
    dim <- dim(y);
    if (is.null(dim)) {
      nbrOfArrays <- length(getDataSet(this));
      return(list(theta=rep(NA, nbrOfArrays),
                  sdTheta=rep(NA, nbrOfArrays),
                  thetaOutliers=rep(NA, nbrOfArrays),
                  phi=c(),
                  sdPhi=c(),
                  phiOutliers=c()
                 )
            );
    }

    if (length(dim) != 2) {
      str(y);
      stop("Argument 'y' must have two dimensions: ",
                                                paste(dim, collapse="x"));
    }

    K <- dim[1];  # Number of probes
    I <- dim[2];  # Number of arrays

    # Too many probes?
    if (K > skipThreshold[1] && I > skipThreshold[2]) {
      warning("Ignoring a unit group when fitting probe-level model, because it has a ridiculously large number of data points: ", paste(dim, collapse="x"), " > ", paste(skipThreshold, collapse="x"));

      naValue <- as.double(NA);
      return(list(theta=rep(naValue, times=I),
                  sdTheta=rep(naValue, times=I),
                  thetaOutliers=rep(naValue, times=I),
                  phi=rep(naValue, times=K),
                  sdPhi=rep(naValue, times=K),
                  phiOutliers=rep(naValue, times=K)
                 )
            );
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Transform data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Add shift
    y <- y + shift;

    # Log-additive model
    y <- log(y, base=2);

    # Look for cells that have NAs in at least one sample?
    w <- NULL;
    nasRemoved <- FALSE;
    if (treatNAsAs == "ignore") {
      hasNAs <- FALSE;
    } else {
      isNA <- is.na(y);
      hasNAs <- any(isNA);
      if (hasNAs) {
        if (treatNAsAs == "weights") {
          badCells <- colAlls(isNA);
          if (any(badCells)) {
            return(list(theta=rep(NA, I),
                        sdTheta=rep(NA, I),
                        thetaOutliers=rep(NA, I),
                        phi=rep(NA, K),
                        sdPhi=rep(NA, K),
                        phiOutliers=rep(NA, K)
                       )
                  );
          }
          w <- matrix(1, nrow=K, ncol=I);
          w[isNA] <- 0;
          y[isNA] <- 0;
        } else if (treatNAsAs == "0") {
          y[isNA] <- 0;
          hasNAs <- FALSE;
        } else if (treatNAsAs == "NA") {
          K0 <- K;  # Number of cells
          okCells <- !rowAnys(isNA);
          # Analyze only valid cells
          y <- y[okCells,,drop=FALSE];
          nasRemoved <- TRUE;
          hasNAs <- FALSE;

          # No valid cells left?
          if (nrow(y) == 0) {
            return(list(theta=rep(NA, I),
                        sdTheta=rep(NA, I),
                        thetaOutliers=rep(NA, I),
                        phi=rep(NA, K0),
                        sdPhi=rep(NA, K0),
                        phiOutliers=rep(NA, K0)
                       )
                  );
          }
        }
      } # if (hasNAs)
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Fit model
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    hasPriors <- !is.null(priors);
    if (hasPriors) {
      paf <- priors$probeAffinities[[1]];
      # Sanity check
      stopifnot(!is.null(paf));
      phi <- paf$phi;
      # Sanity check
      stopifnot(!is.null(phi));
      # Sanity check
      stopifnot(length(phi) == K);

      sdPhi <- paf$sdPhi;

      # Log-additive model
      alpha <- log(phi, base=2);

      if (!is.null(w)) {
        fit <- wrlm(y, alpha, psiCode, psiK, w);
      } else {
        fit <- rlm(y, alpha, psiCode, psiK);
      }
    } else {
      # Use median polish for large probesets (that doesn't have NAs)?
      if (K > medianPolishThreshold[1] && I > medianPolishThreshold[2] && !hasNAs) {
        mp <- medpolish(y, trace.iter=FALSE);
        fit <- list(
          Estimates = c(mp$overall+mp$col, mp$row),
          StdErrors = rep(0, length(c(mp$row, mp$col)))
        );
      } else {
        # Fit model using preprocessCore/affyPLM code
        if (!is.null(w)) {
          fit <- wrlm(y, psiCode, psiK, w);
        } else {
          fit <- rlm(y, psiCode, psiK);
        }
      }
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extract parameters
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    est <- fit$Estimates;
    se <- fit$StdErrors;

    # Chip effects
    beta <- est[1:I];
    # On the intensity scale
    theta <- 2^beta;

    # Probe affinities
    if (!hasPriors) {
      alpha <- est[(I+1):length(est)];
      alpha[length(alpha)] <- -sum(alpha[1:(length(alpha)-1)]);
      # On the intensity scale
      phi <- 2^alpha;
    }

    # The RMA model is fitted with constraint sum(alpha) = 0, that is,
    # such that prod(phi) = 1.

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fit function must return: theta, sdTheta, thetaOutliers,
    # phi, sdPhi, phiOutliers.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (is.null(se)) {
      # For affyPLM v1.10.0 (2006-09-26) or older.
      sdTheta <- rep(1, times=I);
      if (!hasPriors) {
        sdPhi <- rep(1, times=K);
      }
    } else {
      # For affyPLM v1.11.6 (2006-11-01) or newer.
      sdTheta <- 2^(se[1:I]);
      if (!hasPriors) {
        sdPhi <- 2^(se[(I+1):length(se)]);
      }
    }

    # Handle NAs?
    if (nasRemoved) {
      if (treatNAsAs == "NA") {
        naValue <- as.double(NA);
        phi0 <- rep(naValue, times=K0);
        phi0[okCells] <- phi;
        phi <- phi0;

        sdPhi0 <- rep(naValue, times=K0);
        sdPhi0[okCells] <- sdPhi;
        sdPhi <- sdPhi0;

        K <- K0;
      }
    }

    thetaOutliers <- rep(FALSE, times=I);
    phiOutliers <- rep(FALSE, times=K);

    # Return data on the intensity scale
    list(theta=theta, sdTheta=sdTheta, thetaOutliers=thetaOutliers,
         phi=phi, sdPhi=sdPhi, phiOutliers=phiOutliers);
  } # rmaModelAffyPlm()
  attr(rmaModelAffyPlm, "name") <- "rmaModelAffyPlm";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # rmaModelOligo().
  # Author: Henrik Bengtsson, WEHI, 2006-12-11.
  # Requires: oligo() by Benilto Carvalho et al.
  # Why: To fully imitate CRLMM in oligo.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (flavor == "oligo") {
    # First, try to see if package is available;
    pkg <- "oligo";
    .require <- require
    .require(pkg, character.only=TRUE) || throw("Package not loaded: ", pkg)
    pkgDesc <- packageDescription(pkg);
    ver <- pkgDesc$Version;
    verbose && cat(verbose, pkg, " version: ", ver);
    if (compareVersion(ver, "1.7.19") >= 0) {
      # HB 2009-05-09:
      # The API of the native function rma_c_complete_copy()
      # has been updated yet again, but the good thing,
      # there is now a basicRMA() wrapper function.
      # Lookup oligo::basicRMA() once; '::' is expensive
      oligo_basicRMA <- oligo::basicRMA;
      fitRma <- function(y, unitNames, nbrOfUnits, ...) {
        oligo_basicRMA(y, pnVec=unitNames, background=FALSE,
                                            normalize=FALSE, verbose=FALSE);
      } # fitRma()
    } else {
      throw("Non-supported version (< 1.7.19) of 'oligo' detected. Please update the 'oligo' package: ", ver);
    }
  }

  rmaModelOligo <- function(y, priors=NULL, ...) {
    if (!is.null(priors)) {
      throw("NOT IMPLEMENTED: Internal rmaModelOligo() does not support prior parameters.");
    }

    # Add shift
    y <- y + shift;

    # Assert right dimensions of 'y'.
    dim <- dim(y);
    if (length(dim) != 2) {
      str(y);
      stop("Argument 'y' must have two dimensions: ",
                                                paste(dim, collapse="x"));
    }

    K <- dim[1];  # Number of probes
    I <- dim[2];  # Number of arrays

    # Too many probes?
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

    # make factor variables for chip and probe
    unitNames <- rep("X", K);  # dummy probe names
    nbrOfUnits <- as.integer(1); # Only one unit group is fitted

    # Each call to fitRma() outputs "Calculating Expression".
    fit <- fitRma(y, unitNames, nbrOfUnits);

    # Extract probe affinities and chip estimates
    est <- fit[1,,drop=TRUE];  # Only one unit

    # Chip effects
    beta <- est[1:I];

    # Probe affinities
    alpha <- rep(0, K);  # Not returned by fitRma()!

    # Estimates on the intensity scale
    theta <- 2^beta;
    phi <- 2^alpha;

    # The RMA model is fitted with constraint sum(alpha) = 0, that is,
    # such that prod(phi) = 1.

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fit function must return: theta, sdTheta, thetaOutliers,
    # phi, sdPhi, phiOutliers.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    sdTheta <- rep(1, I);
    thetaOutliers <- rep(FALSE, I);
    sdPhi <- rep(1, K);
    phiOutliers <- rep(FALSE, K);

    # Return data on the intensity scale
    list(theta=theta, sdTheta=sdTheta, thetaOutliers=thetaOutliers,
         phi=phi, sdPhi=sdPhi, phiOutliers=phiOutliers);
  } # rmaModelOligo()
  attr(rmaModelOligo, "name") <- "rmaModelOligo";



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Getting the PLM fit function");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the flavor of fitting algorithm for the RMA PLM
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Selecting fit function depending on 'flavor'");
  verbose && cat(verbose, "Flavor: ", flavor);

  priors <- getListOfPriors(this);
  hasPriors <- !is.null(priors);
  verbose && cat(verbose, "Has priors: ", hasPriors);
  if (hasPriors) {
    # Sanity check
    stopifnot(is.element("probeAffinities", names(priors)));
  }

  # Shift signals?
  shift <- this$shift;
  if (is.null(shift))
    shift <- 0;
  verbose && cat(verbose, "Amount of shift: ", shift);

  # Handle non-positive signals?
  treatNAsAs <- this$treatNAsAs;
  if (is.null(treatNAsAs))
    treatNAsAs <- "ignore";
  verbose && cat(verbose, "treatNAsAs: ", treatNAsAs);

  if (flavor == "affyPLM") {
    fcnList <- RmaPlm$getRlmFitFunctions(withPriors=hasPriors, verbose=less(verbose));
    verbose && str(verbose, fcnList);
    # To please R CMD check
    rlm <- wrlm <- NULL; rm(list=c("rlm", "wrlm"));
    attachLocally(fcnList);
    rmaModel <- rmaModelAffyPlm;
  } else if (flavor == "oligo") {
    requireNamespace("oligo") || throw("Package not loaded: oligo")
    rmaModel <- rmaModelOligo;
  } else {
    throw("Cannot get fit function for RMA PLM. Unknown flavor: ", flavor);
  }
  verbose && str(verbose, rmaModel);
  verbose && exit(verbose);

  # Test that it works and is available.
  verbose && enter(verbose, "Validating the fit function on some dummy data");
  ok <- FALSE;
  tryCatch({
    y <- matrix(1:12+0.1, ncol=3);
    if (hasPriors) {
      verbose && cat(verbose, "With prior probe affinities");
      phi <- c(0.7630639, 0.9068485, 1.1208795, 1.2892744);
      sdPhi <- c(1.101706, 1.091220, 1.091220, 1.094554);
      priors <- list(probeAffinities=list("foo"=list(phi=phi, sdPhi=sdPhi)));
      res <- rmaModel(y, priors=priors);
    } else {
      res <- rmaModel(y);
    }
    ok <- TRUE;
  }, error = function(ex) {
    print(ex);
  })
  if (!ok) {
    throw("The fit function for requested RMA PLM flavor failed: ", flavor);
  }
  verbose && exit(verbose);


  verbose && exit(verbose);

  rmaModel;
}, private=TRUE) # getFitUnitGroupFunction()


setMethodS3("getCalculateResidualsFunction", "RmaPlm", function(static, ...) {
  function(y, yhat) {
    y/yhat;
  }
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2012-09-05
# o Now getRlmFitFunctions() of RmaPlm no longer generates a warning
#   in R CMD check about external .Call():s.
# 2012-09-01
# o BUG FIX: getFitUnitGroupFunction() for RmaPlm would throw exception
#   'Error in pkgDesc$Version : $ operator is invalid for atomic vectors'
#   if RmaPlm(..., flavor="oligo") and 'oligo' is not installed.  Now it
#   gives a more informative error.
# 2012-08-21
# o CLEANUP: Dropped support for obsolete RmaPlm(..., flavor="affyPLMold").
# o CLEANUP: Now RmaPlm(..., flavor="oligo") no longer supports
#   oligo (< 1.7.19), which can be considered obsolete version by now.
# o CLEANUP: getRlmFitFunctions() no longer falls back to affyPLM, which
#   was really only needed back in 2007.
# 2012-01-14
# o Updated getFitUnitGroupFunction() for RmaPlm to support prior
#   probe-affinity parameters.
# o Added argument withPriors=FALSE to getRlmFitFunctions() for RmaPlm.
# o ROBUSTNESS: Now the fit functions of RmaPlm and MbeiPlm give an
#   error whenever trying to use prior parameters, which are yet
#   not supported.
# 2009-05-09
# o Updated getFitUnitGroupFunction() of RmaPlm to work with the new
#   oligo v1.7.19 as well, which luckily gut oligo::basicRMA().
# 2008-12-08
# o Now the fit function returned by getFitUnitFunction() must be able
#   to handle prior parameters as well.
# 2008-12-04
# o BUG FIX: flavor='oligo' stopped working awhile ago. Changed the call to
#   'fitRma()' to the native C code as given by
#   'selectMethod("rma", "FeatureSet")'.
# 2008-07-13
# o BUG FIX: plm$treatNAsAs=="NA" returned an incorrect number of probe
#   affinities whenever missing values were exluded.
# 2008-02-12
# o Added mechanism to avoid fitting unit groups with ridiculously many cells.
# 2007-10-06
# o Now the asterisk tag ('*') is no longer assigned in the constructor,
#   but in getTags().
# 2007-09-18
# o BUG FIX: Due to a migration of code from affyPLM to preprocessCore,
#   the fit function returned by getFitFunction() would not work with
#   affyPLM >= 1.13.9.  Now getFitFunction() adopts to the version of
#   affyPLM installed.
# o Updated getFitFunction() with verbose output.
# 2007-09-16
# o Renamed the variables such that index I is for samples and K is for
#   probes, as in the paper.
# 2007-09-15
# o Now the RmaPlm fit function detects cases where a probe get weight zero
#   for all arrays. In such (rare) cases, parameter estimates equals NAs.
# 2007-04-15
# o Added first support for weights in fit function.  This requires
#   affyPLM v1.11.14.  Thanks Ben Bolstad for this.
# 2007-02-14
# o Added getCalculateResidualsFunction().
# 2006-12-12
# o Confirmed that flavor="oligo" replicates the estimates of the oligo
#   package perfectly.
# o Added argument 'flavor' to constructor to make it possible to specify
#   what model fitting algorithm to be used.  If other than the default
#   flavor, a tag is added to indicate what flavor is used.
# 2006-11-02
# o Added SE estimates in RmaPlm from Ben's new code. Works with
#   affyPLM v1.11.6 or newer. /KS
# 2006-09-26
# o Added code to use either of the two RMA fit functions.
# o Incorporated Ken Simpson's fit function for RMA as an alternative.
# 2006-09-11
# o The fit function now returns all required fields.
# 2006-08-25
# o Created from the corresponding Li & Wong model.
############################################################################
