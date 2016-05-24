
  C1C2toAB <- function(xy, ...) {
    dxy <- colDiffs(xy[,1:2]);
#print(dxy);
    # (l,b) = Length and slope of line
    l <- sqrt(dxy[,1]^2 + dxy[,2]^2);
    k <- (dxy[,2]/dxy[,1]);
    b <- atan(1/k);
#print(list(l=l, b=b, k=k, dxy=dxy));


    # Regions weights
    counts <- xy[,"dhNbrOfLoci", drop=TRUE];
    w <- sqrt(counts);
    w <- w / sum(w, na.rm=TRUE);

    # Change-point weights
    w <- w[1:(length(w)-1)] + w[2:length(w)];

    res <- as.matrix(cbind(l=l, b=b, w=w, k=k));
#    attr(res, "dxy") <- dxy;
    res;
  } # C1C2toAB()

  ABtoC1C2 <- function(AB, C1C2, ...) {
    C1C2 <- C1C2[,c("C1","C2"), drop=FALSE];
    C1C2 <- as.matrix(C1C2);
    dxy <- colDiffs(C1C2);

    #  b = atan(dC1C2[,1]/dC1C2[,2]);
    #  tan(b) = dC1/dC2;
    #  dC2*tan(b) = dC1;
    l <- AB[,"l"];
    b <- AB[,"b"];
    k <- AB[,"k"];
    dxy2 <- sign(dxy);
    k2 <- tan(b);
    dxy2[,2] <- dxy2[,1]/k2;
    l2 <- sqrt(dxy2[,1]^2+dxy2[,2]^2);
    s <- l/l2;
    dxy2 <- s*dxy2;

    for (cc in 1:2) {
      delta <- dxy2[,cc];
      idxs <- which(is.finite(delta));
      C1C2[idxs+1,cc] <- C1C2[idxs,cc] + delta[idxs];
    }

    C1C2;
  } # ABtoC1C2()

###########################################################################/**
# @set "class=PairedPSCBS"
# @RdocMethod normalizeBAFsByRegions
#
# @title "Normalizes allele B fractions (BAFs) based on region-based PSCN estimates"
#
# \description{
#  @get "title" as given by the PSCBS segmentation method.
# }
#
# @synopsis
#
# \arguments{
#   \item{fit}{A PairedPSCBS fit object as returned by
#     @see "PSCBS::segmentByPairedPSCBS".}
#   \item{by}{A @character string specifying if the normalization function
#     should be estimated based on TumorBoost normalized or non-normalized
#     tumor allele B fractions (BAFs).}
#   \item{...}{Additional arguments passed
#     @see "aroma.cn::normalizeMirroredBAFsByRegions".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a PairedPSCBS fit object where the region-level
#   decrease-in-heterozygosity (DH) means have been normalized,
#   as well as the locus-specific tumor allele B fractions.
# }
#
# \details{
#   Note that his normalization method depends on the segmentation
#   results. Hence, it recommended \emph{not} to resegment the
#   normalized signals returned by this, because such a segmentation
#   will be highly dependent on the initial segmentation round.
# }
#
# @examples "../incl/normalizeBAFsByRegions.PairedPSCBS.Rex"
#
# @author "HB, PN"
#
# \seealso{
#   Internally @see "aroma.cn::normalizeMirroredBAFsByRegions" is used.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("normalizeBAFsByRegions", "PairedPSCBS", function(fit, by=c("betaTN", "betaT"), ..., force=FALSE, cache=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'by':
  by <- match.arg(by);

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'cache':
  cache <- Arguments$getLogical(cache);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Normalizes region-level mirrored allele B fractions (mBAFs)");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  key <- list(method="normalizeBAFsByRegions", class=class(fit)[1],
    data=as.data.frame(fit),
    version="2011-11-02"
  );
  dirs <- c("aroma.cn", "ortho");
  if (!force) {
    res <- loadCache(key=key, dirs=dirs);
    if (!is.null(res)) {
      verbose && cat(verbose, "Cached results found.");
      verbose && exit(verbose);
      return(res);
    }
  }


  data <- getLocusData(fit);
  stopifnot(!is.null(data));

  segs <- getSegments(fit, splitters=TRUE);
  stopifnot(!is.null(segs));

  chromosomes <- getChromosomes(fit);
  nbrOfChromosomes <- length(chromosomes);
  verbose && cat(verbose, "Number of chromosomes: ", nbrOfChromosomes);
  verbose && print(verbose, chromosomes);

  nbrOfSegments <- nrow(segs);
  verbose && cat(verbose, "Number of segments: ", nbrOfSegments);

  # Get mean estimators
  estList <- getMeanEstimators(fit, "dh");
  avgDH <- estList$dh;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chromosome <- data$chromosome;
  x <- data$x;
  betaT <- data$betaT;
  betaTN <- data$betaTN;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Calculate region-level mBAFs for homozygous SNPs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating region-level mBAFs for homozygous SNPs");

  # Calculate mBAFs for all loci
  beta <- data[[by]];
  rho <- 2*abs(beta - 1/2);
  # Not needed anymore
  beta <- NULL;

  # Identify homozygous SNPs
  muN <- data$muN;
  isHom <- (muN == 0 | muN == 1);
  # Not needed anymore
  muN <- NULL;

  # Drop all homozygous SNPs already here.
  # TO DO

  # Allocate
  naValue <- as.double(NA);
  X <- matrix(naValue, nrow=nbrOfSegments, ncol=3);
  for (kk in seq_len(nbrOfSegments)) {
    chrKK <- as.numeric(segs[kk,"chromosome"]);
    xRange <- as.numeric(segs[kk,c("dhStart", "dhEnd")]);
    tcn <- segs[kk,"tcnMean"];
    dh <- segs[kk,"dhMean"];
    # Identify all homozygous SNPs in the region
    keep <- (chromosome == chrKK & xRange[1] <= x & x <= xRange[2] & isHom);
    keep <- which(keep);
    mBAFhom <- avgDH(rho[keep], na.rm=TRUE);
    X[kk,] <- c(dh, mBAFhom, tcn);
  } # for (kk ...)

  # Not needed anymore
  rho <- isHom <- NULL;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize region-level DHs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing region-level mBAFs");
  XN <- normalizeMirroredBAFsByRegions(data=X, ..., verbose=verbose);

  # Update DH segmentation means
  rhoN <- XN[,1,drop=TRUE];
  segs[,"dhMean"] <- rhoN;
  # Not needed anymore
  rhoN <- NULL;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Normalize locus-level data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Normalizing locus-level data accordingly");
  modelFit <- attr(XN, "modelFit");
  scale <- modelFit$scale;
  # Not needed anymore
  modelFit <- XN <- NULL;

  # Expand region-level scale factors to locus-level scale factors
  naValue <- as.double(NA);
  scales <- rep(naValue, times=length(data$betaT));
  for (kk in seq_len(nbrOfSegments)) {
    chrKK <- as.numeric(segs[kk,"chromosome"]);
    xRange <- as.numeric(segs[kk,c("dhStart", "dhEnd")]);
    # Identify all SNPs in the region
    keep <- (chromosome == chrKK & xRange[1] <= x & x <= xRange[2]);
    keep <- which(keep);
    scales[keep] <- scale[kk];
  } # for (kk ...)

  # Update tumor allele B fractions
  for (ff in c("betaT", "betaTN")) {
    beta <- data[[ff]];
    beta <- beta - 1/2;
    beta <- scales * beta;
    beta <- beta + 1/2;
    data[[ff]] <- beta;
  }
  # Not needed anymore
  beta <- scales <- NULL;
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fitN <- fit;
  fitN$data <- data;
  fitN$output <- segs;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    saveCache(key=key, dirs=dirs, fitN);
  }

  verbose && exit(verbose);

  fitN;
}) # normalizeBAFsByRegions()



##############################################################################
# HISTORY
# 2013-01-17 [HB]
# o Updated normalizeBAFsByRegions() for PairedPSCBS to recognize when
#   other mean-level estimators than the sample mean have been used.
# 2011-10-16 [HB]
# o Now using getLocusData(fit) and getSegments(fit) where applicable.
# 2011-07-10 [HB]
# o Updated code to work with the new column names in PSCBS v0.11.0.
# 2010-10-10 [HB]
# o Added memoization to normalizeBAFsByRegions().
# 2010-09-26 [HB]
# o Now normalizeBAFsByRegions() for PairedPSCBS handles multiple chromosomes.
# 2010-09-19 [HB+PN]
# o Added orthogonalizeC1C2() for PairedPSCBS.
# 2010-09-15 [HB]
# o Added Rdocs for callCopyNeutralRegions().
# 2010-09-09 [HB]
# o Added callCopyNeutralRegions() for PairedPSCBS.
# 2010-09-08 [HB]
# o Added subsetBySegments() for PairedPSCBS.
# o Added Rdocs with an example.
# o Added normalizeBAFsByRegions() for PairedPCSBS.
# o Created.
##############################################################################
