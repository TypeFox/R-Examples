###########################################################################/**
# @RdocClass NmfPlm
#
# @title "The NmfPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the NMF model of [REF].
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "aroma.affymetrix::ProbeLevelModel".}
#   \item{maxIter}{The maximum number of iteration in the NMF step.}
#   \item{maxIterRlm}{A positive @integer specifying the maximum number of
#     iterations used in rlm.}
#   \item{refs}{An index @vector (@integer or @logical) specifying the
#     reference samples.  If @NULL, all samples are used as a reference.}
#   \item{flavor}{(Internal/developmental only)
#      A @character string specifying which algorithm to use.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \seealso{
#   Internally, for each SNP the NMF model is fitted using the
#   @see "fitSnpNmf" function.
# }
#
# \references{
#  [1] @include "../incl/OrtizM_etal_2010.Rd" \cr
# }
#
# @author
#*/###########################################################################
setConstructorS3("NmfPlm", function(..., maxIter=10L, maxIterRlm=20L, refs=NULL, flavor=c("v4", "v3", "v2", "v1")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'maxIter':
  maxIter <- Arguments$getInteger(maxIter, range=c(1L,999L));

  # Argument 'maxIterRlm':
  maxIterRlm <- Arguments$getInteger(maxIterRlm, range=c(1L,999L));

  # Argument 'flavor':
  flavor <- match.arg(flavor);


  extend(ProbeLevelModel(...), "NmfPlm",
    .maxIter = maxIter,
    .maxIterRlm = maxIterRlm,
    .refs = refs,
    .flavor = flavor
  )
})



setMethodS3("getAsteriskTags", "NmfPlm", function(this, collapse=NULL, ...) {
  # Returns 'PLM[,<shift>]'

  tags <- NextMethod("getAsteriskTags", collapse=NULL);
  tags[1] <- "NMF";

  flavor <- this$.flavor;
  if (!is.null(flavor) && flavor != "v4") {
    tags <- c(tags, flavor);
  }

  # Collapse
  tags <- paste(tags, collapse=collapse);

  tags;
})



setMethodS3("getProbeAffinityFile", "NmfPlm", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get the probe affinities (and create files etc)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paf <- NextMethod("getProbeAffinityFile");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Update the encode and decode functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  setEncodeFunction(paf, function(groupData, ...) {
    phi <- .subset2(groupData, "phiA");
    stdvs <- .subset2(groupData, "phiB");
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
      phiA=intensities,
      phiB=stdvs,
      phiOutliers=outliers
    );
  })
  paf;
}, private=TRUE)


setMethodS3("getFitUnitFunction", "NmfPlm", function(this,...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get algorithm parameters
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the amount of shift to be used.
  shift <- this$shift;
  if (is.null(shift))
    shift <- 0;

  # Algorithm version
  flavor <- this$.flavor;
  if (is.null(flavor)) {
    flavor <- "v4";
  }

  # Maximum number of iterations to fit.
  maxIter <- this$.maxIter;
  if (is.null(maxIter)) {
    maxIter <- 10L;
  }

  # Maximum number of iterations to fit rlm.
  maxIterRlm <- this$.maxIterRlm;
  if (is.null(maxIterRlm)) {
    maxIterRlm <- 20L;
  }

  # Reference samples.
  refs <- this$.refs;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Select the 'NMF' function to use
  # (When adding a new version, add it here; not below)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  nmfFcn <- NULL;
  if (flavor == "v1") {
    nmfFcn <- get("NmfPinv", mode="function");
  } else if (flavor == "v2") {
    nmfFcn <- get("NmfFast", mode="function");
  } else if (flavor == "v3") {
    nmfFcn <- get("NmfPinvBeta", mode="function");
  } else if (flavor == "v4") {
    nmfFcn <- get("fitSnpNmf", mode="function");
  } else {
    throw("Unknown flavor: ", flavor);
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup the fit function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitUnits <- function(y, prevpaf=NULL,...){
    SNPunit <- y;
    #save(SNPunit, file = "SNPunit")
    groupNames <- names(SNPunit);

    if (length(SNPunit) > 1L) {
      SNPdata <- rbind(SNPunit[[1L]]$intensities, SNPunit[[2L]]$intensities);
      SNPdata <- SNPdata + shift;
      indNonPos <- (SNPdata <= 0);
      SNPdata[indNonPos] <- 0.0001;
      nbrGroups <- 2L;
    } else {
      SNPdata <- SNPunit[[1L]]$intensities;
      nbrGroups <- 1L;
    }

    if (nbrGroups == 2L) {
      NMFdata <- nmfFcn(SNPdata, maxIter=maxIter, maxIterRlm=maxIterRlm, refs=refs);

      W <- NMFdata[[1L]];
      H <- NMFdata[[2L]];

      I <- dim(H)[2L];
      K <- dim(W)[1L]/2;

      # prepare returned data
      # allele A
      theta1 <- H[1L,];
      sdTheta <- rep(1, times=I);
      thetaOutliers <- logical(I);
      phi1 <- W[1:K,1L];
      sdPhi1 <- W[(K+1):(2*K),1L];
      phiOutliers <- logical(K);

      # allele B
      theta2 <- H[2L,];
      phi2 <- W[1:K,2L];
      sdPhi2 <- W[(K+1):(2*K),2L];

      # fitted unit
      fitUU <- list(
        A = list(theta=theta1, sdTheta=sdTheta, thetaOutliers=thetaOutliers, phiA=phi1, phiB=sdPhi1, phiOutliers=phiOutliers),
        B = list(theta=theta2, sdTheta=sdTheta, thetaOutliers=thetaOutliers, phiA=phi2, phiB=sdPhi2, phiOutliers=phiOutliers)
      );
    } else {
      I <- dim(SNPdata)[2L];
      K <- dim(SNPdata)[1L];
      theta1 <-  rep(1, times=I);
      sdTheta <- rep(1, times=I);
      thetaOutliers <- logical(I);

      phi1 <- double(K);
      sdPhi1 <- double(K);
      phiOutliers <- logical(K);

      fitUU <- list(list(theta=theta1, sdTheta=sdTheta, thetaOutliers=thetaOutliers, phiA=phi1, phiB=sdPhi1, phiOutliers=phiOutliers));
    }

    names(fitUU) <- groupNames;
    fitUU;
  } # getFitUnitFunction()


  fitUnits;
}, private=TRUE)


############################################################################
# HISTORY:
# 2010-09-28 [HB]
# o Now argument 'refs' defaults to NULL (not 0), which means all samples.
# 2010-06-04 [MO]
# o Added refs as argument.
# 2010-05-18 [MO]
# o Added maxIterRlm as argument.
# 2010-05-17 [HB]
# o Now a flavor tag is added to NmfPlm:s only if flavor != "v4" (default).
# 2009-11-18 [MO]
# o Removed internal save() in getFitUnitFunction() of NmfPlm.
# 2009-03-24 [HB]
# o Added Rdoc comments.
# 2009-01-28 [HB]
# o Made getFitUnitFunction() slightly faster. Cleaned up code. Added
#   support for 'flavor' to specify which NMF fitting function to use.
# 2008-12-08 [HB]
# o Tidied up code.
# o Updated to make use of new ProbeLevelModel.R.
############################################################################
