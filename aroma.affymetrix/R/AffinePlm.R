###########################################################################/**
# @RdocClass AffinePlm
# \encoding{latin1}
#
# @title "The AffinePlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents affine model in Bengtsson \& Hossjer (2006).
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "ProbeLevelModel".}
#   \item{background}{If @TRUE, background is estimate for each unit group,
#     otherwise not. That is, if @FALSE, a \emph{linear} (proportional)
#     model without offset is fitted, resulting in very similar results as
#     obtained by the @see "MbeiPlm".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   For a single unit group, the affine model is:
#
#    \deqn{y_{ik} = a + \theta_i \phi_k + \varepsilon_{ik}}
#
#   where \eqn{a} is an offset common to all probe signals,
#   \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I},
#   and \eqn{\phi_k} are the probe affinities for probes \eqn{k=1,...,K}.
#   The \eqn{\varepsilon_{ik}} are zero-mean noise with equal variance.
#   The model is constrained such that \eqn{\prod_k \phi_k = 1}.
#
#   Note that with the additional constraint \eqn{a=0} (see arguments above),
#   the above model is very similar to @see "MbeiPlm".  The differences in
#   parameter estimates is due to difference is assumptions about the
#   error structure, which in turn affects how the model is estimated.
# }
#
# @author "HB"
#
# \references{
#   Bengtsson \& Hossjer (2006). \cr
# }
#*/###########################################################################
setConstructorS3("AffinePlm", function(..., background=TRUE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load required packages
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...);
  if (length(args) > 0 && !is.null(args[[1]])) {
    # Early error, iff package is missing
    requireNamespace("aroma.light") || throw("Package not loaded: aroma.light");
  }

  this <- extend(ProbeLevelModel(...), "AffinePlm",
    background = background
  );
  validate(this);
  this;
})


setMethodS3("getAsteriskTags", "AffinePlm", function(this, collapse=NULL, ...) {
  # Returns 'PLM[,<shift>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "AFF";

  # Add class specific parameter tags
  if (!this$background)
    tags <- c(tags, "lin");

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)


setMethodS3("getProbeAffinityFile", "AffinePlm", function(this, ...) {
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
})



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
# @author
#
# \seealso{
#   @see "aroma.light::calibrateMultiscan".
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitUnitGroupFunction", "AffinePlm", function(this, ...) {
  standardize <- this$standardize;
  center <- this$background;
  shift <- this$shift;
  if (is.null(shift))
    shift <- 0;

  affineFit <- function(y, ...) {
    # Add shift
    y <- y + shift;

    # NOTE: If center=FALSE => constraint a=0 /HB 2006-09-11
    y <- t(y);
    f <- .calibrateMultiscan(y, center=center, project=TRUE);
    theta <- as.vector(f);
    phi <- as.vector(attr(f, "modelFit")$b);

    I <- length(theta);
    K <- length(phi);

    # Rescale such that prod(phi) = 1?
    if (standardize) {
      c <- prod(phi)^(1/K);
      phi <- phi/c;
      theta <- theta*c;
    }

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
  } # affineFit()

  affineFit;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2007-12-06
# o The tag for the affine PLM is now 'AFF' (before it was APLM).
# 2007-10-06
# o Now the asterisk tag ('*') is no longer assigned in the constructor,
#   but in getTags().
# 2007-09-16
# o Renamed the variables such that index I is for samples and K is for
#   probes, as in the paper.
# 2007-03-29
# o Changed tag 'linear' to 'lin'.
# 2006-09-11
# o Added argument 'background' to fit background or not.
# 2006-08-28
# o Created from the Li & Wong model.
############################################################################
