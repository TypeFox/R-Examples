# Estimates factor scores and parameters for PLS path models
sempls <- function(model, ...){
  UseMethod("sempls", model)
}

sempls.plsm <-
function(model, data, maxit=20, tol=1e-7, scaled=TRUE, sum1=FALSE, wscheme="centroid", pairwise=FALSE,
         method=c("pearson", "kendall", "spearman"),
         convCrit=c("relative", "square"), verbose=TRUE, ...){
  method <- match.arg(method)
  convCrit <- match.arg(convCrit)
  result <- list(coefficients=NULL, path_coefficients=NULL,
                 outer_loadings=NULL ,cross_loadings=NULL,
                 total_effects=NULL,inner_weights=NULL, outer_weights=NULL,
                 blocks=NULL, factor_scores=NULL, data=NULL, scaled=scaled,
                 model=model, weighting_scheme=NULL, weights_evolution=NULL,
                 sum1=sum1, pairwise=pairwise, method=method, iterations=NULL,
                 convCrit=convCrit, verbose=verbose, tolerance=tol, maxit=maxit, N=NULL,
                 incomplete=NULL, Hanafi=NULL)
  class(result) <- "sempls"

  # checking the data
  data <- data[, model$manifest]
  N <- nrow(data)
  missings <- which(complete.cases(data)==FALSE)
  if(length(missings)==0 & verbose){
    cat("All", N ,"observations are valid.\n")
    if(pairwise){
      pairwise <- FALSE
      cat("Argument 'pairwise' is reset to FALSE.\n")
    }
  }
  else if(length(missings)!=0 & !pairwise & verbose){
    # Just keeping the observations, that are complete.
    data <- na.omit(data[, model$manifest])
    cat("Data rows:", paste(missings, collapse=", "),
        "\nare not taken into acount, due to missings in the manifest variables.\n",
        "Total number of complete cases:", N-length(missings), "\n")
  }
  else if(verbose){
     cat("Data rows", paste(missings, collapse=", "),
         " contain missing values.\n",
         "Total number of complete cases:", N-length(missings), "\n")
  }
  ## check the variances of the data
  if(!all(apply(data, 2, sd, na.rm=TRUE) != 0)){
     stop("The MVs: ",
          paste(colnames(data)[which(apply(data, 2, sd)==0)], collapse=", "),
          "\n  have standard deviation equal to 0.\n",
          "  Recheck model!\n")
  }

  ## scale data?
  # Note: scale() changes class(data) to 'matrix'
  if(scaled) data <- scale(data)

  ## compute PLS approximation of LV scores

  #############################################
  # step 1: Initialisation
  stp1 <- step1(model, data, sum1=sum1, pairwise, method)
  factor_scores <- stp1$latent
  if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale"),
                  each=length(model$manifest))
      Wold <- stp1$outerW / sdYs
  }
  else Wold <- stp1$outerW
  weights_evolution <- reshape(as.data.frame(Wold),
                               v.names="weights",
                               ids=rownames(Wold),
                               idvar="MVs",
                               times=colnames(Wold),
                               timevar="LVs",
                               varying=list(colnames(Wold)),
                               direction="long")
  ## fixme: find something more efficient than 'cbind'
  weights_evolution <- cbind(weights_evolution, iteration=0)
  Hanafi <- cbind(f=sum(abs(cor(factor_scores)) * model$D),
                  g=sum(cor(factor_scores)^2 * model$D),
                  iteration=0)


  #############################################
  # Select the function according to the weighting scheme
  if(wscheme %in% c("A", "centroid")) {
    innerWe <- centroid
    result$weighting_scheme <- "centroid"
  }
  else if(wscheme %in% c("B", "factorial")) {
    innerWe <- factorial
    result$weighting_scheme <- "factorial"
  }
  else if(wscheme %in% c("C", "pw", "pathWeighting")) {
    innerWe <- pathWeighting
    result$weighting_scheme <- "path weighting"
  }
  else {stop("The argument E can only take the values 'A', 'B' or 'C'.\n See ?sempls")}

  converged <- c()
  i <- c()
  Wnew <- c()
  innerWeights <- c()
  eval(plsLoop)

  ## print
  if(converged & verbose){
      cat(paste("Converged after ", (i-1), " iterations.\n",
                "Tolerance: ", tol ,"\n", sep=""))
      if (wscheme %in% c("A", "centroid")) cat("Scheme: centroid\n")
      if (wscheme %in% c("B", "factorial")) cat("Scheme: factorial\n")
      if (wscheme %in% c("C", "pw", "pathWeighting")) cat("Scheme: path weighting\n")
  }
  else if(!converged){
      stop("Result did not converge after ", result$maxit, " iterations.\n",
           "\nIncrease 'maxit' and rerun.", sep="")
  }

  weights_evolution <- weights_evolution[weights_evolution!=0,]
  weights_evolution$LVs <- factor(weights_evolution$LVs,  levels=model$latent)
  # create result list
  ifelse(pairwise, use <- "pairwise.complete.obs", use <- "everything")
  result$path_coefficients <- pathCoeff(model=model, factor_scores, method, pairwise)
  result$cross_loadings <- cor(data, factor_scores, use, method)
  result$outer_loadings <- result$cross_loadings
  result$outer_loadings[Wnew==0] <- 0
  result$total_effects <- totalEffects(result$path_coefficients)
  result$inner_weights <- innerWeights
  result$outer_weights <- Wnew
  result$weights_evolution <- weights_evolution
  result$Hanafi <- Hanafi
  result$factor_scores <- factor_scores
  result$data <- data
  result$N <- N
  result$incomplete <- missings
  result$iterations <- (i-1)
  result$coefficients <- coefficients(result)
  return(result)
}

plsLoop <- expression({
  #######################################################################
  # Iterate over step 2 to step 5
  i <- 1
  converged <- FALSE
  while(!converged){

    #############################################
    # step 2
    innerWeights <- innerWe(model, fscores=factor_scores, pairwise, method)
    factor_scores <- step2(Latent=factor_scores, innerWeights, model, pairwise)

    #############################################
    # step 3
    Wnew <-  outerApprx2(Latent=factor_scores, data, model,
                        sum1=sum1, pairwise, method)

    #############################################
    # step 4
    factor_scores <- step4(data, outerW=Wnew, model, pairwise)
    if(!sum1){
      # to ensure: w'Sw=1
      sdYs <- rep(attr(factor_scores, "scaled:scale"),
                  each=length(model$manifest))
      Wnew <- Wnew / sdYs
    }
    weights_evolution_tmp <- reshape(as.data.frame(Wnew),
                                     v.names="weights",
                                     ids=rownames(Wnew),
                                     idvar="MVs",
                                     times=colnames(Wnew),
                                     timevar="LVs",
                                     varying=list(colnames(Wnew)),
                                     direction="long")
    weights_evolution_tmp <- cbind(weights_evolution_tmp, iteration=i)
    weights_evolution <- rbind(weights_evolution, weights_evolution_tmp)
    Hanafi_tmp <- cbind(f=sum(abs(cor(factor_scores)) * model$D),
                         g=sum(cor(factor_scores)^2 * model$D),
                         iteration=i)
    Hanafi <- rbind(Hanafi, Hanafi_tmp)

    #############################################
    # step 5
    st5 <- step5(Wold, Wnew, tol, converged, convCrit)
    Wold <- st5$Wold
    converged <- st5$converged

    #############################################


    if(i == maxit && !converged){
      # 'try-error' especially for resempls.R
      class(result) <- c(class(result), "try-error")
      i <- i+1
      break
    }

    i <- i+1
  }
})

# print method
print.sempls <- function(x, digits=2, ...){
  print(x$coefficients, digits=digits, ...)
  invisible(x)
}
