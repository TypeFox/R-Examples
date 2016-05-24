
# the plan here is to shuffle the ranefs back into the way a merMod object
# stores them so that a simple X * beta + Z * theta op does the trick


fitted.sim.merMod <- function(object, regression,...){
  if (missing(regression) || is.null(regression)) stop("fitted for sim.mer requires original merPred object as well.");
  if (!inherits(regression, "merMod")) stop("regression argument for fitted on sim.mer does not inherit from class 'merMod'");
  sims <- object;
  numSimulations <- dim(sims@fixef)[1];
  devcomp <- getME(regression, "devcomp");
  dims <- devcomp$dims;

  numRanef  <- dims[["q"]];
  numLevels <- dims[["reTrms"]];

  simulatedRanef <- matrix(0, numRanef, numSimulations);

  index <- 0;
  for (i in 1:length(sims@ranef)) {
    levelSims <- sims@ranef[[i]];
    numCoefficientsPerLevel <- dim(levelSims)[2];
    numGroupsPerLevel <- dim(levelSims)[3];
    for (j in 1:numCoefficientsPerLevel) {
      ranefRange <- index + 1:numGroupsPerLevel;
      index <- index + numGroupsPerLevel;

      simulatedRanef[ranefRange,] <- t(levelSims[,j,]);
    }
  }

  X <- getME(regression, "X");
  Zt <- getME(regression, "Zt");

  linearPredictor <- as.matrix(tcrossprod(X, sims@fixef) + crossprod(Zt, simulatedRanef)) +
    matrix(getME(regression, "offset"), dims[["n"]], numSimulations);

  if (dims[["GLMM"]] == 0L){
    return(linearPredictor)
  }else{
      return(regression@resp$family$linkinv(linearPredictor))
  }
};
