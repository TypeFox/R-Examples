###########################################################################/**
# @RdocClass AvgPlm
#
# @title "The AvgPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents a PLM where the probe intensities are averaged
#  assuming identical probe affinities.
#  For instance, one may assume that replicated probes with identical
#  sequences have the same probe affinities, cf. the GenomeWideSNP\_6
#  chip type.
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
#   For a single unit group, the averaging PLM of K probes is:
#
#    \deqn{y_{ik} = \theta_i + \varepsilon_{ik}}
#
#   where \eqn{\theta_i} are the chip effects for arrays \eqn{i=1,...,I}.
#   The \eqn{\varepsilon_{ik}} are zero-mean noise with equal variance.
# }
#
# \section{Different flavors of model fitting}{
#   The above model can be fitted in two ways, either robustly or
#   non-robustly.
#   Use argument \code{flavor="mean"} to fit the model non-robustly, i.e.
#
#    \deqn{\hat{\theta}_{i} = 1/K \sum_k y_{ik}}.
#
#   Use argument \code{flavor="median"} to fit the model robustly, i.e.
#
#    \deqn{\hat{\theta}_{i} = median_k y_{ik}}.
#
#   Missing values are always excluded.
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("AvgPlm", function(..., flavor=c("median", "mean")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'flavor':
  flavor <- match.arg(flavor);

  this <- extend(ProbeLevelModel(...), "AvgPlm",
    .flavor = flavor
  );
  validate(this);
  this;
})


setMethodS3("validate", "AvgPlm", function(this, ...) {
  ds <- getDataSet(this);
  if (is.null(ds))
    return(invisible(TRUE));

  if (length(ds) < 1) {
    throw("This ", class(this)[1], " requires at least 1 array: ",
                                                         length(ds));
  }

  invisible(TRUE);
}, protected=TRUE)


setMethodS3("getAsteriskTags", "AvgPlm", function(this, collapse=NULL, ...) {
  # Returns 'PLM[,<shift>]'
  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "AVG";

  # Add class specific parameter tags
  if (this$.flavor != "median")
    tags <- paste(tags, this$.flavor, sep=",");

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
}, protected=TRUE)



setMethodS3("getParameters", "AvgPlm", function(this, ...) {
  params <- NextMethod("getParameters");
  params$flavor <- this$.flavor;
  params;
}, protected=TRUE)



setMethodS3("getProbeAffinityFile", "AvgPlm", function(this, ...) {
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
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFitUnitGroupFunction", "AvgPlm", function(this, ...) {
  # Float precision
#  .Machine$float.eps <- sqrt(.Machine$double.eps);
  floatEps <- sqrt(.Machine$double.eps);

  # Shift signals?
  shift <- this$shift;
  if (is.null(shift))
    shift <- 0;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # avgModel()
  # Author: Henrik Bengtsson, UC Berkeley.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  avgPlmModel <- function(y, ...){
    # Assert right dimensions of 'y'.

    # If input data are dimensionless, return NAs.
    if (is.null(dim(y))) {
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

    if (length(dim(y)) != 2) {
      str(y);
      stop("Argument 'y' must have two dimensions: ",
                                                paste(dim(y), collapse="x"));
    }

    # Add shift
    y <- y + shift;

    I <- ncol(y);  # Number of arrays
    K <- nrow(y);  # Number of probes

    # Fit model
    if (K == 1) {
      theta <- y;
      sdTheta <- rep(0,I); # floatEps;
    } else {
      y <- t(y);
      if (flavor == "median") {
        theta <- rowMedians(y, na.rm=TRUE);
        sdTheta <- rowMads(y, center=theta, na.rm=TRUE);
      } else if (flavor == "mean") {
        theta <- .rowMeans(y, na.rm=TRUE);
        sdTheta <- rowSds(y, mean=theta, na.rm=TRUE);
      }
    }

    # Should we store std deviations or std errors?!? /HB 2007-09-08
    sdTheta <- sdTheta/sqrt(K);

    # Hmm..., when searching for "units todo", we use
    # (sdTheta <= 0).
    sdTheta <- sdTheta + floatEps;

    # Probe affinities are all identical (==ones)
    phi <- rep(1, K);
    sdPhi <- rep(1, K);  # Default, we not estimated (should we store NAs?!?)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fit function must return: theta, sdTheta, thetaOutliers,
    # phi, sdPhi, phiOutliers.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    thetaOutliers <- rep(FALSE, I);
    phiOutliers <- rep(FALSE, K);

    # Return data on the intensity scale
    list(theta=theta, sdTheta=sdTheta, thetaOutliers=thetaOutliers,
         phi=phi, sdPhi=sdPhi, phiOutliers=phiOutliers);
  } # avgPlmModel()
  attr(avgPlmModel, "name") <- "avgPlmModel";



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the flavor of fitting algorithm for the averaging PLM
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  flavor <- this$.flavor;
  if (flavor == "median") {
    fitFcn <- avgPlmModel;
  } else if (flavor == "mean") {
    fitFcn <- avgPlmModel;
  } else {
    throw("Cannot get fit function for AvgPlm. Unknown flavor: ", flavor);
  }

  # Test that it works and is available.
  ok <- FALSE;
  tryCatch({
    fitFcn(matrix(1:6+0.1, ncol=3));
    ok <- TRUE;
  }, error = function(ex) {
    print(ex);
  })

  if (!ok) {
    throw("The fit function for requested AvgPlm flavor failed: ", flavor);
  }

  fitFcn;
}, private=TRUE)


setMethodS3("getCalculateResidualsFunction", "AvgPlm", function(static, ...) {
  function(y, yhat) {
    y - yhat;
  }
}, static=TRUE, protected=TRUE)



############################################################################
# HISTORY:
# 2008-12-31
# o Added an ad hoc validate() to AvgPlm to overriden the default test
#   in MultiarrayUnitModel (yes, because it is still not a true single
#   array implementation).
# 2007-10-06
# o Added getAsteriskTag() to AvgPlm.
# 2007-09-16
# o Renamed the variables such that index I is for samples and K is for
#   probes, as in the paper.
# 2007-09-12
# o WORKAROUND: If there is only one probe, the fit function will return
#   theta=y:s, and sdTheta=0:s.  However, when searching for units to do,
#   we test (sdTheta <= 0).  The workaround is to store the smallest
#   float available instead of zero.
# 2007-09-08
# o Created from RmaPlm.R.
############################################################################
