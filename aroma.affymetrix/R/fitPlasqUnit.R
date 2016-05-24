# Gets PLASQ probe types from an AffymetrixCdfFile.
setMethodS3("readPlasqTypes", "AffymetrixCdfFile", function(this, ..., verbose=FALSE) {
  pathname <- getPathname(this);
  cdf <- .readCdf(pathname, ..., readUnitDirection=TRUE, readGroupDirection=TRUE, readXY=FALSE, readIndexpos=FALSE, readIsPm=TRUE, readAtoms=FALSE, readUnitAtomNumbers=FALSE, readGroupAtomNumbers=FALSE);
  cdf <- getPlasqTypes(cdf);
  cdf;
}, private=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Redefine gaussian to allow for exponential link.
plasqGaussian <- function(link="exponential") {
  linktemp <- substitute(link);
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp);
    if (linktemp == "link")
      linktemp <- eval(link);
  }

  if (linktemp == "exponential") {
    stats <- list(
      linkfun = function(mu) exp(mu),
      linkinv = function(eta) log(eta),
      mu.eta = function(eta) 1/eta,
      valideta = function(eta) all(eta > 0)
    );
  }

 structure(list(
   family = "gaussian",
   link = linktemp,
   linkfun = stats$linkfun,
   linkinv = stats$linkinv,
   variance = function(mu) {
     rep.int(1, length(mu))
   },
   dev.resids = function(y, mu, wt) {
     wt * ((y - mu)^2)
   },
   aic = function(y, n, mu, wt, dev) {
     sum(wt) * (log(dev/sum(wt) * 2 * pi) + 1) + 2
   },
   mu.eta = stats$mu.eta,
   initialize = expression({
     n <- rep.int(1, nobs)
     mustart <- y
   }),
   validmu = function(mu) TRUE
 ), class = "family");
} # plasqGaussian()


# Adopted and optimized from EMSNP() in PLASQ500K.
setMethodS3("fitPlasqUnit", "matrix", function(ly, ptype, maxIter=1000, acc=0.1, ..., glmFamily=plasqGaussian("exponential"), .digits=NULL) {


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Main
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dim <- dim(ly);
  nbrOfProbes <- dim[1];
  nbrOfSamples <- dim[2];
  sampleNames <- colnames(ly);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Indices of different probe types:
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # 'ptype' coding:
  #   0=MMoBR,  1=MMoBF,  2=MMcBR,  3=MMcBF,  4=MMoAR,  5=MMoAF,
  #   6=MMcAR,  7=MMcAF,  8=PMoBR,  9=PMoBF, 10=PMcBR, 11=PMcBF,
  #  12=PMoAR, 13=PMoAF, 14=PMcAR, 15=PMcAF.
  isPm <- (ptype %in% 8:15);
  isA <- (ptype %in% c(4:7,12:15));
  isFwd <- (ptype %in% c(1,3,5,7,9,11,13,15));
  isCentered <- (ptype %in% c(2,3,6,7,10,11,14,15));

#  # Create 'ptype'
#  ptype0 <- ptype;
#  ptype <- rep(NA, nbrOfProbes);
#  ptype[!isPm & !isA & !isFwd & !isCentered] <-  0;
#  ptype[!isPm & !isA &  isFwd & !isCentered] <-  1;
#  ptype[!isPm & !isA & !isFwd &  isCentered] <-  2;
#  ptype[!isPm & !isA &  isFwd &  isCentered] <-  3;
#  ptype[!isPm &  isA & !isFwd & !isCentered] <-  4;
#  ptype[!isPm &  isA &  isFwd & !isCentered] <-  5;
#  ptype[!isPm &  isA & !isFwd &  isCentered] <-  6;
#  ptype[!isPm &  isA &  isFwd &  isCentered] <-  7;
#  ptype[ isPm & !isA & !isFwd & !isCentered] <-  8;
#  ptype[ isPm & !isA &  isFwd & !isCentered] <-  9;
#  ptype[ isPm & !isA & !isFwd &  isCentered] <- 10;
#  ptype[ isPm & !isA &  isFwd &  isCentered] <- 11;
#  ptype[ isPm &  isA & !isFwd & !isCentered] <- 12;
#  ptype[ isPm &  isA &  isFwd & !isCentered] <- 13;
#  ptype[ isPm &  isA & !isFwd &  isCentered] <- 14;
#  ptype[ isPm &  isA &  isFwd &  isCentered] <- 15;
#  stopifnot(identical(ptype, ptype0));

  pmA <- which(isPm &  isA);
  pmB <- which(isPm & !isA);
  mmC <- which(!isPm & isCentered);
  mmAo <- which(!isPm &  isA & !isCentered);
  mmBo <- which(!isPm & !isA & !isCentered);

  isFwd <- rep(isFwd, nbrOfSamples);
  fIdxs <- which(isFwd);
  rIdxs <- which(!isFwd);
  nbrOfFwd <- length(fIdxs);
  nbrOfRev <- length(rIdxs);
  hasFwd <- (nbrOfFwd > 0);
  hasRev <- (nbrOfRev > 0);

  ### Matrix that gives parameter indices for each ptype (row-1)
  # A 16x4 matrix with values in [1,13].
  paramIndMat <- cbind(
    rep(c(6,1),8),
    c(11,11,rep(c(9,4),5),7,2,7,2),
    c(10,5,10,5,11,11,10,5,8,3,8,3,10,5,10,5),
    rep(c(13,12),8)
  );
  colnames(paramIndMat) <- c("gamma", "alpha", "beta", "sigma");
  paramIndMat2 <- paramIndMat[ptype+1,];
  gammaIdxs <- paramIndMat2[,1];
  alphaIdxs <- paramIndMat2[,2];
  betaIdxs  <- paramIndMat2[,3];
  sigmaIdxs <- paramIndMat2[,4];
  # Not needed anymore
  paramIndMat <- paramIndMat2 <- NULL;

  iis <- seq_len(nbrOfSamples);
  offsets <- nbrOfProbes*(iis-1);

  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  # E_0-STEP (Initialization)
  # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
  # Samples: i = 1,2,...,N
  # SNP: j0
  # Copy-number levels: l = 0,1,2 (==BB,AB,AA)
  # Z_il = I(C_Ai = l) = indicator for C_Ai == l.
  # z_il = E[Z_il]

  # Note: here we use z_li, not z_il.
  naValue <- as.double(NA);
  zs <- matrix(naValue, nrow=3, ncol=nbrOfSamples);
  ok <- !is.na(ly[1,]);
  for (ii in iis[ok]) {
    # H0: muA > muB
    # H1: muA <= muB
    lyA <- ly[pmA, ii];
    lyB <- ly[pmB, ii];

    # P value
    p <- t.test(lyA, lyB, "greater")$p.value;

    if (p > 0.5) {
      zs[,ii] <- c((3*p-1)/2, 1-p, (1-p)/2);
    } else {
      zs[,ii] <- c(p/2, p, 1-3*p/2);
    }
  } # for (ii in ...)
#  verfy(zs, key=list("zs"));


  # Vectorize
  ly <- as.vector(ly);

  # Intial parameter estimates
  betasF <- c(0,0,0,0,0);
  betasR <- c(0,0,0,0,0);
  ps <- c(1,1,0);

  # Allocate variables
  mu <- vector("list", 3);
  ms <- double(3);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Expectation-Maximation algorithm
  #
  # Treating {Z_il}_il as latent variables.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  converged <- FALSE;
  for (rr in seq_len(maxIter)) {
    zs.old <- zs;
    betasF.old <- betasF;
    betasR.old <- betasR;
    ps.old <- ps;

    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # M-STEP
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # Appendix p17:
    #   \hat{p}_{j0} = \frac{1}{N} \sum_{i=1}^N z_{ij_0l}
    # Proportion of samples with genotypes (BB, AB, AA)
    ps <- rowSums(zs) / nbrOfSamples;

#    verfy(ps, key=list(iter=rr, "ps"));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Setup design matrix
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # CA =          Z_i1 + 2*Z_i2; i=1,2,...,N
    # CB = 2*Z_i0 + Z_i1         ; i=1,2,...,N
    CA <-            zs[2,] + 2*zs[3,];
    CB <- 2*zs[1,] + zs[2,];

    # Columns: (PM_A, PM_B, ..., ...)
    X <- matrix(0, nrow=nbrOfProbes*nbrOfSamples, ncol=4);
    for (ii in iis) {
      if (!is.na(zs[1,ii])) {
        offset <- offsets[ii];
        X[offset + pmA              , 1] <- CA[ii];
        X[offset + pmB              , 2] <- CB[ii];
        X[offset + c(pmB, mmC, mmAo), 3] <- CA[ii];
        X[offset + c(mmC, pmA, mmBo), 4] <- CB[ii];
      }
    } # for (ii in ...)

#    verfy(X, key=list(iter=rr, "X"));

#    verfy(fIdxs, key=list(iter=rr, "Find"));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # A fully homozygote SNP?
    if (min(zs[1,], na.rm=TRUE) > .9*.7 || min(zs[3,], na.rm=TRUE) > .9*.7) {
      X[,1] <- X[,1]+X[,2];
      X[,2] <- X[,3]+X[,4];
      X <- X[,1:2];
#    verfy(X, key=list(iter=rr, "X", "sub"));

      if (hasFwd) {
        fit <- glm(ly[fIdxs] ~ X[fIdxs,], family=glmFamily,
                                                       na.action=na.omit);
        betasF <- coef(fit)[c(1,2,2,3,3)];
        sigF <- sqrt(sum(residuals(fit, "deviance")^2)/(nbrOfFwd-3));
      } else {
        betasF <- rep(NA, 5);
        sigF <- NA;
      }

      if (hasRev) {
        fit <- glm(ly[rIdxs] ~ X[rIdxs,], family=glmFamily,
                                                         na.action=na.omit);
        betasR <- coef(fit)[c(1,2,2,3,3)];
        sigR <- sqrt(sum(residuals(fit, "deviance")^2)/(nbrOfRev-3));
      } else {
        betasR <- rep(NA, 5);
        sigR <- NA;
      }
    } else {
      if (hasFwd) {
        fit <- glm(ly[fIdxs] ~ X[fIdxs,], family=glmFamily,
                                                         na.action=na.omit);
        betasF <- coef(fit)
        sigF <- sqrt(sum(residuals(fit, "response")^2)/(nbrOfFwd-5));
      } else {
        betasF <- rep(NA, 5);
        sigF <- NA;
      }

      if (hasRev) {
        fit <- glm(ly[rIdxs] ~ X[rIdxs,], family=glmFamily,
                                                        na.action=na.omit);
        betasR <- coef(fit);
        sigR <- sqrt(sum(residuals(fit, "response")^2)/(nbrOfRev-5));
      } else {
        betasR <- rep(NA,5);
        sigR <- NA;
      }
    }

#params<-c(betasF, betasR, 0, sigF, sigR)
#attributes(params) <- NULL;
#    verfy(params, key=list(iter=rr, "params"));


    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    # E-step
    # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Extract (gamma, alpha, beta, sigma) for each probe
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # 'params' is a vector of length 13.
    # params = (
    #   gammaF, alphaF0, betaF0, alphaF1, betaF1,
    #   gammaR, alphaR0, betaR0, alphaR1, betaR1,
    #   0, sigmaF, sigmaR
    # )
    params <- c(betasF, betasR, 0, sigF, sigR);
    gamma <- params[gammaIdxs];
    alpha <- params[alphaIdxs];
    beta  <- params[betaIdxs];
    sigma <- params[sigmaIdxs];

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Calculate the means for each of the genotypes
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    mu[[1]] <- log(gamma +           2*beta);
    mu[[2]] <- log(gamma +   alpha +   beta);
    mu[[3]] <- log(gamma + 2*alpha         );

##     mu00 <- list(mu0, mu1, mu2);
##     mu00 <- lapply(mu00, FUN=unname);
##
##     # The PLASQ implementation of the above two steps:
##     # For each probe, select the four parameters.
##     # Input: 16x4 matrix. Output: 4x16 matrix
##     pMat <- apply(paramIndMat2, MARGIN=1, FUN=function(idxs) {
##       params[idxs];
##     });
##     pMat <- t(pMat);
##
##     # Appendix p16:
##     # Pre-calculate some constants
##     mu0 <- log((pMat %*% c(1,0,2,0))[,1]);
##     mu1 <- log((pMat %*% c(1,1,1,0))[,1]);
##     mu2 <- log((pMat %*% c(1,2,0,0))[,1]);
##     # Not needed anymore
##     pMat <- params <- NULL; # Not needed anymore
##     mu11 <- list(mu0, mu1, mu2);
##
##     stopifnot(identical(mu00, mu11));


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # For each sample, find the expected genotypes, i.e. z = E[Z]
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    for (ii in iis) {
      if (!is.na(zs[1,ii])) {
        offset <- offsets[ii];
        yy <- ly[offset + 1:nbrOfProbes];

        # Appendix p16: f_0(.), f_1(.), f_2(.).
        for (kk in 1:3)
          ms[kk] <- prod(dnorm(yy, mean=mu[[kk]], sd=sigma), na.rm=TRUE);

        # Appendix p18: z_ij = p_l * f_l(y_i) / f_{Y_i}(y_i)
        # Appendix p16: f_Y(.) = sum_{l=0}^2 p_l f_l(.)
        zs[,ii] <- ps * ms / (ps %*% ms);
      }
    }

    # Check for convergance
    diffs <- max(abs(c(
      (zs-zs.old), (ps-ps.old), (betasF-betasF.old), (betasR-betasR.old)
    )), na.rm=TRUE);
    converged <- (diffs <= acc);
    if (converged)
      break;
  } # for (rr in ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Model parameters
  names(betasF) <- c("gamma_F", "alpha_F0", "beta_F0", "alpha_F1", "beta_F1");
  names(betasR) <- c("gamma_R", "alpha_R0", "beta_R0", "alpha_R1", "beta_R1");

  # Genotypes
  rownames(zs) <- c("P[BB]", "P[AB]", "P[AA]");
  colnames(zs) <- sampleNames;

  # CA =          Z_i1 + 2*Z_i2; i=1,2,...,N
  # CB = 2*Z_i0 + Z_i1         ; i=1,2,...,N
  CA <-            zs[2,] + 2*zs[3,];
  CB <- 2*zs[1,] + zs[2,];
  CN <- matrix(c(CA,CB), nrow=2, byrow=TRUE);
  rownames(CN) <- c("A", "B");

  if (!is.null(.digits)) {
    zs <- round(zs, digits=.digits);
    CN <- round(CN, digits=.digits);
  }

  list(
    coeffsF=betasF,
    coeffsR=betasR,
    sigs=c(sigF, sigR),
    z=zs,
    CN=CN,
    converged=converged,
    nbrOfIterations=rr
  );
}, private=TRUE) # fitPlasqUnit()


############################################################################
# HISTORY:
# 2008-07-20
# o Updated the following methods to preallocate matrixes with the correct
#   data type to avoid coercing later: fitPlasqUnit().
# 2007-01-11
# o Replaces argument 'mat' (on intensity scale) with 'ly' (on log scale).
# 2007-01-01
# o Verified that this implementation gives identical results to
#   PLASQ500K::EMSNP().
# 2006-12-31
# o Created fitPlasqUnit() to fit a single CEL unit.  Got code to get the
#   PLASQ probe types from a CDF too.  This is major requirement for
#   implementing PLASQ in aroma.affymetrix. TO DO: Verify correctness.
# 2006-12-29
# o Created.
############################################################################
