###########################################################################/**
# @RdocClass MbeiPlm
#
# @title "The MbeiPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the \emph{model-based expression indexes} (MBEI)
#  multiplicative model in Li \& Wong (2001).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   For a single unit group, the multiplicative model is:
#
#    \deqn{y_{ik} = \theta_i \phi_k + \varepsilon_{ik}}
#
#   where \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I},
#   and \eqn{\phi_k} are the probe affinities for probes \eqn{k=1,...,K}.
#   The \eqn{\varepsilon_{ik}} are zero-mean noise with equal variance.
#   To make to parameters identifiable, the constraint
#   \eqn{\prod_k \phi_k = 1} is added.
# }
#
# @author "HB"
#
# \seealso{
#   Internally @see "affy::fit.li.wong" is used.
# }
#
# \references{
#   Li, C. and Wong, W.H. (2001), Genome Biology 2, 1-11.\cr
#   Li, C. and Wong, W.H. (2001), Proc. Natl. Acad. Sci USA 98, 31-36.\cr
# }
#*/###########################################################################
setConstructorS3("MbeiPlm", function(...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...);
  if (length(args) > 0 && !is.null(args[[1]])) {
    # Early error, iff package not installed.
    requireNamespace("affy") || throw("Package not loaded: affy");
  }

  this <- extend(ProbeLevelModel(...), "MbeiPlm");
  validate(this);
  this;
})


setMethodS3("getAsteriskTags", "MbeiPlm", function(this, collapse=NULL, ...) {
  # Returns 'PLM[,<shift>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "MBEI";

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getProbeAffinityFile", "MbeiPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the probe affinities (and create files etc)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paf <- NextMethod("getProbeAffinityFile");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the encode and decode functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setEncodeFunction(paf, function(groupData, ...) {
    # Rename some fields so that we support the structure of this class,
    # but also output from affy::fit.li.wong().
    names <- names(groupData);
    # Is it an affy:fit.li.wong() structure?
    if ("sdPhi" %in% names) {
      names <- sub("iter", "nbrOfIterations", names);
      names <- sub("convergence1", "converged", names);
      names <- sub("convergence2", "convergedOutliers", names);
      names(groupData) <- names;
    }

    # Encode outliers as the sign of 'pixels'; -1 = TRUE, +1 = FALSE
    pixels <- sign(0.5 - as.integer(groupData$phiOutliers));

    # Encode the number of iterations as the absolute value of the 1st pixel.
    pixels[1] <- pixels[1]*groupData$nbrOfIterations;

    # Encode convergence1/2 as bits in the 2nd pixel.
    # Note: For this to work, there must be at least two 'pixels' in the
    # unit group.  However, this is not true for some strands on the 500K
    # chip, which might have a single probe quartet, e.g. SNP_A-1781633 on
    # the Nsp chip.  Thus, this is a problem if the MBEI model is fitted
    # strand by strand.  On the other hand, the PLM based on chip effects
    # and probe affinities breaks down if there is only one probe, so
    # talking about probe affinities does not make sense in the first
    # place.  Thus, the probe affinity for this single probe will be set
    # to exactly one (on the intensity scale) by the construct that the
    # probe affinities should multiply to one.  There is no standard
    # deviation estimate or convergence results for such single-probe
    # models. /HB 2006-12-18
    npixels <- length(pixels);
    if (npixels > 1) {
      pixels[2] <- pixels[2] *
              (1 + 2*groupData$converged + 4*groupData$convergedOutliers);
    }

    list(intensities=groupData$phi, stdvs=groupData$sdPhi, pixels=pixels);
  })

  setDecodeFunction(paf,  function(groupData, ...) {
    pixels <- groupData$pixels;

    # Outliers are encoded by the sign of 'pixels'.
    outliers <- as.logical(1-sign(pixels));

    # Number of iterations is encoded as the absolute value of the 1st pixel.
    nbrOfIterations <- as.integer(abs(pixels[1])+0.5);

    npixels <- length(pixels);
    if (npixels > 1) {
      # convergence & convergenceOutliers are encoded as bits in the 2nd pixel.
      t <- pixels[2] %/% 2;
      converged <- as.logical(t %% 2 == 1);  t <- t %/% 2;
      convergedOutliers <- as.logical(t %% 2 == 1);
    } else {
      # See comments on the encode function above. /HB 2006-12-18
      converged <- TRUE;
      convergedOutliers <- TRUE;
    }

    list(
      phi=groupData$intensities,
      sdPhi=groupData$stdvs,
      phiOutliers=outliers,
      nbrOfIterations=nbrOfIterations,
      converged=converged,
      convergedOutliers=convergedOutliers
    );
  })

  paf;
}, private=TRUE)



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
# }
#
# \value{
#  Returns a @function.
# }
#
# \seealso{
#   @see "affy::fit.li.wong".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitUnitGroupFunction", "MbeiPlm", function(this, ...) {
  # To help 'R CMD check' to locate fit.li.wong().
  requireNamespace("affy") || throw("Package not loaded: affy");
  fit.li.wong <- affy::fit.li.wong


  standardize <- this$standardize;
  shift <- this$shift;
  if (is.null(shift))
    shift <- 0;

  liWong <- function(y, priors=NULL, ...) {
    if (!is.null(priors)) {
      throw("NOT IMPLEMENTED: Internal liWong() does not support prior parameters.");
    }

    # Add shift
    y <- y + shift;

    # Enough of probes?
    if (nrow(y) > 1) {
      fit <- fit.li.wong(t(y));

      # A fit function must return: theta, sdTheta, thetaOutliers,
      # phi, sdPhi, phiOutliers.
      names <- names(fit);
      idxs <- match(c("sigma.theta", "theta.outliers", "sigma.phi",
                                                     "phi.outliers"), names);
      names[idxs] <- c("sdTheta", "thetaOutliers", "sdPhi", "phiOutliers");
      names(fit) <- names;

      # Rescale such that prod(phi) = 1?
      if (standardize) {
        phi <- fit$phi;
        theta <- fit$theta;
        K <- length(phi);
        c <- prod(phi)^(1/K);
        phi <- phi/c;
        theta <- theta*c;
        fit$phi <- phi;
        fit$theta <- theta;
      }
    } else {
      # For the case where there is only a single probe in the unit group
      # let the chip effect be the probe signal and the probe affinity one.
      # This is indeed what fit.li.wong() returns, but we don't want their
      # warning about it:
      fit <- list(theta=as.vector(y), sdTheta=NA, thetaOutliers=NA,
                  phi=1, sdPhi=NA, phiOutliers=NA,
                  sigma.eps=NA, single.outliers=NA,
                  convergence1=NA, convergence2=NA,
                  iter=1, delta=NA);
    }

    fit;
  } # liWong()

  liWong;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2012-03-23
# o Now getFitUnitGroupFunction() for MbeiPlm helps 'R CMD check' to
#   locate fit.li.wong() by explicitly requiring 'affy'.
# 2012-01-14
# o ROBUSTNESS: Now the fit functions of RmaPlm and MbeiPlm give an
#   error whenever trying to use prior parameters, which are yet
#   not supported.
# 2006-12-18
# o Now the fit function of the MBEI model treats single-probe unit groups
#   specially; the affy::fit.li.wong() handled it already before, but
#   generated a warning for each call.
# o BUG FIX: The encoding/decoding of probe affinities assumed at least two
#   probes per unit group, but this is not true for all chip types, e.g.
#   the 500K SNP chips.
# 2006-09-10
# o Updated getFitFunction() to return required fields.
# 2006-08-24
# o Added Rdoc comments.
# 2006-08-23
# o Added getProbeAffinities() and the corrsponding cached fields.
# o Now fit() does not re-read data just updated.
# 2006-08-19
# o After all the bug fixes in updateCel() I think this function finally
#   works correctly.
# 2006-08-17
# o Created.
############################################################################
