# simulations of sigma, fixef, and ranef drawn from a posterior
# under a flat prior and conditioned on estimate of ranef covar
setMethod("sim", signature(object = "merMod"),
          function(object, n.sims=100)
{
  applyLeftFactor <- function(decomp, rhs) {
    c(as.vector(decomp$ul %*% rhs[ranefRange] + decomp$ur %*% rhs[fixefRange]),
      as.vector(decomp$lr %*% rhs[fixefRange]));
  }

  # information is conditional on hyperparameters
  # information is of [ranef, fixef]
  getInverseInformationLeftFactor <- function(regression) {
    Lz  <- getME(regression, "L");
    Rzx <- getME(regression, "RZX");
    Rx  <- getME(regression, "RX");

    # upper left, lower right, and lower left blocks of left-factor
    # of inverse
    solveFunc <- getMethod("solve", signature(a = "CHMfactor", b = "diagonalMatrix"));
    Rz.inv  <- t(solveFunc(Lz, Diagonal(Lz@Dim[1]), "L"));
    Rx.inv  <- solve(Rx);
    Rzx.inv <- -Rz.inv %*% Rzx %*% Rx.inv;

    # this is me figuring some stuff out. new lmer doesn't permute Zt apparently
    #
    #Lz.tmp <- as(Lz, "sparseMatrix");
    #P.chol <- as(Lz@perm + 1, "pMatrix");
    #Zt <- getME(regression, "Zt");
    #W <- Diagonal(numObs, regression@resp$sqrtXwt);
    ## P.ranef <- getRanefPerm(regression);
    #Lambdat <- getME(regression, "Lambdat") # t(P.ranef) %*% getME(regression, "Lambdat") %*% P.ranef;
    #A <- Lambdat %*% Zt;
    #C <- A %*% W;
    #L.hyp <- Cholesky(tcrossprod(P.chol %*% C), Imult = 1, LDL = FALSE, perm = FALSE);
    #L.hyp@perm <- Lz@perm;
    #L.hyp@type[1] <- 2L;
    #browser();

    #P.ranef <- getRanefPerm(model);
    #Lambda <- P.ranef %*% getRanefChol(model) %*% t(P.ranef);
    Lambda <- t(getME(regression, "Lambda"));

    return(list(ul = Lambda %*% Rz.inv, ur = Lambda %*% Rzx.inv, lr = Rx.inv));
  }

  # assumes p(sigma^2) propto sigma^-2
  sampleCommonScale <- function(ignored) {
    return(sqrt(1 / rgamma(1, 0.5 * numDoF,
                           0.5 * devcomp$cmp[["pwrss"]])));
  }

  regression <- object;
  devcomp <- getME(regression, "devcomp");
  dims <- devcomp$dims;

  if (dims[["NLMM"]] != 0L) stop("sim not yet implemented for nlmms");

  numObs    <- dims[["n"]];
  numRanef  <- dims[["q"]];
  numFixef  <- dims[["p"]];
  numLevels <- dims[["reTrms"]];

  isLinearMixedModel <- dims[["GLMM"]] == 0L && dims[["NLMM"]] == 0L;
  numEffects <- numRanef + numFixef;
  numDoF <- numObs - numFixef;

  # pertain to simulations that we do all as a single vector
  ranefRange <- 1:numRanef;
  fixefRange <- numRanef + 1:numFixef;

  # stuff used to rearrange ranef into usable form
  groupsPerUniqueFactor <- lapply(regression@flist, levels);
  factorPerLevel <- attr(regression@flist, "assign");

  coefficientNamesPerLevel <- regression@cnms;
  numCoefficientsPerLevel <- as.numeric(sapply(coefficientNamesPerLevel, length));
  numGroupsPerLevel <- as.numeric(sapply(groupsPerUniqueFactor[factorPerLevel], length));
  numRanefsPerLevel <- numCoefficientsPerLevel * numGroupsPerLevel;
  ranefLevelMap <- rep.int(seq_along(numRanefsPerLevel), numRanefsPerLevel);

  # storage for sims
  simulatedSD    <- if (isLinearMixedModel) { rep(NA, n.sims); } else { NA };
  simulatedRanef <- vector("list", numLevels);
  names(simulatedRanef) <- names(regression@cnms);
  for (i in 1:numLevels) {
    simulatedRanef[[i]] <- array(NA, c(n.sims, numGroupsPerLevel[i], numCoefficientsPerLevel[i]),
                                 list(NULL, groupsPerUniqueFactor[[factorPerLevel[i]]], coefficientNamesPerLevel[[i]]));
  }

  simulatedFixef <- matrix(NA, n.sims, numFixef,
                           dimnames = list(NULL, names(fixef(regression))));


  # "b" are the rotated random effects, i.e. what ranef() returns in
  # a rearranged format.
  effectsMean <- c(getME(regression, "b")@x, getME(regression, "beta"));
  effectsCovLeftFactor <- getInverseInformationLeftFactor(regression);

  for (i in 1:n.sims) {
    if (isLinearMixedModel) {
      simulatedSD[i] <- sampleCommonScale(regression);
      sphericalEffects <- rnorm(numEffects, 0, simulatedSD[i]);
    } else {
      sphericalEffects <- rnorm(numEffects);
    }
    simulatedEffects <- applyLeftFactor(effectsCovLeftFactor, sphericalEffects) + effectsMean;

    simulatedFixef[i,] <- simulatedEffects[fixefRange];

    rawRanef <- simulatedEffects[ranefRange];
    simulatedRanefPerLevel <- split(rawRanef, ranefLevelMap);
    for (k in 1:numLevels) {
      simulatedRanef[[k]][i,,] <- matrix(simulatedRanefPerLevel[[k]], ncol = numCoefficientsPerLevel[k], byrow = TRUE);
    }
  }

  ans <- new("sim.merMod",
             "fixef" = simulatedFixef,
             "ranef" = simulatedRanef,
             "sigma" = simulatedSD);
  return(ans);
});
