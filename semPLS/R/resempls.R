# Refitting the sempls object for the bootsempls method.
resempls <-
function(sempls, data, start=c("ones", "old"), method, ...){
  start <- match.arg(start)
  bootMethod <- method
  sum1 <- sempls$sum1
  tol <- sempls$tolerance
  maxit <- sempls$maxit
  model <- sempls$model
  pairwise <- sempls$pairwise
  method <- sempls$method
  convCrit <- sempls$convCrit
  result <- sempls

  ## scale data?
  # Note: scale() changes class(data) to 'matrix'
  if(sempls$scaled) data <- scale(data)

  #############################################
  # step 1
  # initialization with results from fitted model
  if(start=="old"){
    factor_scores <- sempls$factor_scores
    Wold <- sempls$outer_weights
    i <- sempls$iterations
    index <- sempls$weights_evolution$iteration==i
    weights_evolution <- sempls$weights_evolution[index,]
    Hanafi <- sempls$Hanafi[i+1,] # Values of the last iteration
  }
  else if(start=="ones"){
    # Weights not adding up to 1 (14.08.2009)
    # changed (30.03.2010)
    stp1 <- step1(model, data, sum1=sum1, pairwise, method)
    factor_scores <- stp1$latent
    if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale"),
                  each=length(model$manifest))
      Wold <- stp1$outerW / sdYs
    }
    else Wold <- stp1$outerW
    index <- sempls$weights_evolution$iteration==0
    weights_evolution <- sempls$weights_evolution[index,]
    Hanafi <- sempls$Hanafi[1,] # iteration==0
  }

  #############################################
  # weighting scheme
  innerWe <- eval(parse(text=sub(" w", "W", sempls$weighting_scheme)))


  #############################################
  # Iterate over step 2 to 5
  i <- c()
  eval(plsLoop)

  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  ### bootstrap method ##################################################
  # Construct level changes
  clcIndex <- NULL
  if(bootMethod=="ConstructLevelChanges"){
    result$cross_loadings <- cor(data, factor_scores, use=use, method=method)
    result$outer_loadings <- result$cross_loadings
    result$outer_loadings[Wnew==0] <- 0
    clcIndex <- which(abs(colSums(sempls$outer_loadings - result$outer_loadings)) >
                    abs(colSums(sempls$outer_loadings + result$outer_loadings)))
    # change the signs of the weights
    if(sum1){
      Wnew[, clcIndex] <- -1*apply(as.matrix(Wnew[, clcIndex]), 2, sum1)
    }
    else {Wnew[, clcIndex] <- -1* Wnew[, clcIndex]}
    
    # repeat step4
    factor_scores <- step4(data, outerW=Wnew, model, pairwise)
    if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale"),
                  each=length(model$manifest))
      Wnew <- Wnew / sdYs
    }
  }
  #######################################################################
  result$cross_loadings <- cor(data, factor_scores, use=use, method=method)
  result$outer_loadings <- result$cross_loadings
  result$outer_loadings[Wnew==0] <- 0
  result$path_coefficients <- pathCoeff(model=model, factor_scores, method, pairwise)
  result$total_effects <- totalEffects(result$path_coefficients)
  result$inner_weights <- innerW
  result$outer_weights <- Wnew
  result$weights_evolution <- weights_evolution
  result$factor_scores <- factor_scores
  result$data <- data
  result$iterations <- (i-1)
  result$coefficients <- coefficients(result)
  result$clcIndex <- clcIndex
  return(result)
}
