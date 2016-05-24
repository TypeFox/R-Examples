`TMC` <-
function(num.Iter=50, data.Size=100, make.Data = gen.Data, make.Params = gen.Params,
    model.List, weight.Vector = rep(1, times = length(model.List)), msc.List, 
    fit.Model = fit.Models, stepSize = .05, sumstats = list("Median Rank" = median ), huge=FALSE, 
    var.Frame=data.frame(), par.Sigma=1, data.Sigma=1, barebones=FALSE, allow.Negs = FALSE,
    thresholds = c(1,2,3,5,10), test.Size = 0, scale.Frame = TRUE, use.Ranks = TRUE, ...){

  ptm <- proc.time()

  # Generates a matrix, each row of which contains weights which together form
  # a convex combination of the MSCs in msc.List.  convCombMat2 is similar,
  # but allows for more general linear combinations, including those with negative
  # weights, and the weights do not necessarily sum to unity.
  
  convCombMat <- function(numFactors = 3, stepSize = stepSize, tol = 10^-9, huge = FALSE){
    if(numFactors == 1) stop("Need more than one MSC")
    if((round(1/stepSize))^numFactors > 500000){if(!huge)
      {stop("Matrix will be BIG... set huge = TRUE to go on anyway.")}}
    tol <- 10^-9
    big <- seq(0,1,by=stepSize)
    l = length(big)
    mat <- matrix(nrow=l^(numFactors-1), ncol = numFactors)
    if(numFactors == 2){mat[,1] = big
      mat[,2] = 1 - big
      return(mat)
    }
    for (i in 1:(numFactors-1)){
      mat[,i] <- rep(big, times = l^(i-1), each = l^(numFactors - 1 - i))
    }
    mat <- zapsmall(mat)
    mat[,numFactors] <- 1 - apply(mat[,1:(numFactors-1)], 1, sum)
    zapsmall(mat[mat[,numFactors] >= 0,])
}
 convCombMat2 <- function(numFactors = 3, stepSize = stepSize, tol = 10^-9, huge = FALSE){
    if(numFactors == 1) stop("Need more than one MSC")
    if((round(1/stepSize))^numFactors > 500000){if(!huge)
      {stop("Matrix will be BIG... set huge = TRUE to go on anyway.")}}
    tol <- 10^-9
    big <- seq(-1,1,by=stepSize)
    l = length(big)
    mat <- matrix(nrow=l^(numFactors), ncol = numFactors)
    if(numFactors == 2){mat[,1] = big
      mat[,2] = 1 - big
      return(mat)
    }
    for (i in 1:(numFactors)){
      mat[,i] <- rep(big, times = l^(i-1), each = l^(numFactors - i))
    }
    mat <- zapsmall(mat)
    mat
}
  if(!is.numeric(test.Size) || test.Size != round(test.Size) || test.Size < 0) stop("test.Size must be a nonnegative integer")
  if(!is.list(msc.List))stop("msc.List must be a list")
  if(!is.list(model.List))stop("model.List must be a list")
  if(!is.list(sumstats))stop("sumstats must be a list")
  if(any(weight.Vector < 0) || !is.numeric(weight.Vector)) stop("weight.Vector must be a vector of non-negative numbers")
  if(!is.function(make.Data))stop("make.Data must be a function")
  if(!is.function(make.Params))stop("make.Params must be a function")
  if(!is.function(fit.Model))stop("fit.Model must be a function")
  if(!is.data.frame(var.Frame))stop("var.Frame must be a data frame")
  
  if(barebones) {
    sumstats <- as.list(rep("Internal Function", times = length(thresholds)+1))
    names(sumstats) <- c("Average Rank", paste("P(Rank > ", thresholds, ")", sep=''))
  }
  if(!(deparse(substitute(var.Frame)) %in% search()))attach(var.Frame)
  on.exit(try(detach("var.Frame"), silent = TRUE))
  if(scale.Frame) var.Frame <- data.frame(scale(var.Frame))
  
  # useful later in plotting functions if all elements of msc.List are named
  namevec <- names(msc.List)
  if(is.null(namevec))namevec <- rep("", length(msc.List))
  namevec[namevec == ""] <- sapply(match.call()$msc.List, deparse)[-1][namevec == ""]
  names(msc.List) = namevec
  
  if(("holdout.Mean" %in% names(msc.List) || "holdout.Med" %in% names(msc.List) || "holdout.SS" %in% names(msc.List))
  && test.Size < 1) stop("Cannot use holdout msc functions if test.Size is not positive!")
  
  test.Vector <- numeric()
  errs <- rep(0, num.Iter)
  weight.Vector <- weight.Vector/sum(weight.Vector)
  num.Crits <- length(msc.List)
  ifelse(allow.Negs, crits.Mat <- convCombMat2(num.Crits, stepSize=stepSize, huge=huge),
    crits.Mat <- convCombMat(num.Crits, stepSize=stepSize, huge=huge))
  num.Stats <- ifelse(!barebones, length(sumstats), length(thresholds))
  simulated.Models <- list()
  simulated.Parameters <- list()
  simulated.Data <- list()
  simulation.Attempts <- rep(0,num.Iter)

  # These will ultimately be parts of the returned msc object.
  if(!barebones) ranks.Mat <- matrix(nrow = dim(crits.Mat)[1], ncol = num.Iter)
  Avg.Rank <- rep(0, dim(crits.Mat)[1]); Probs <- rep(0,dim(crits.Mat)[1]*length(thresholds))
  dim(Probs) <- c(dim(crits.Mat)[1], length(thresholds))
  Sum.Stats <- matrix(nrow = dim(crits.Mat)[1], ncol = num.Crits+num.Stats)
  for(i in 1:num.Iter){

  # in case we need a holdout sample
  if(length(var.Frame) != 0 ){
    if(test.Size != 0){train.Inds <- 1:(dim(var.Frame)[1] - test.Size)
         test.Inds <- (dim(var.Frame)[1] - test.Size + 1):dim(var.Frame)[1]
         test.Frame <- data.frame(var.Frame[test.Inds,])
         train.Frame <- data.frame(var.Frame[train.Inds,])
     }
     else {train.Frame <- var.Frame; test.Frame <- NULL}
  }
  
  else {test.Frame = NULL; train.Frame = NULL}
    
    success.Sim = FALSE
    # Sometimes time-series simulation screws up
    while(success.Sim == FALSE){
      rand <- runif(1)
      
      # Pick a random model from model.List using weight.Vector
      model.Index <- which.max(cumsum(weight.Vector)>rand)
      true.Model <- model.List[[model.Index]]
      
      try(pars <- make.Params(true.Model, par.Sigma=par.Sigma, var.Frame = var.Frame, ...),
      silent=TRUE)

      try({dat <- make.Data(fmla=true.Model, pars=pars,
      data.Size=data.Size, var.Frame = var.Frame, data.Sigma=data.Sigma,
      train.Frame = train.Frame, test.Size = test.Size, ...);
      data.Vector <- dat$data.Vector;
      test.Vector <- dat$test.Vector}, silent=TRUE )
      simulation.Attempts[i] <- simulation.Attempts[i]+1
      if(exists("data.Vector") && is.numeric(data.Vector) && sum(!is.na(data.Vector))>0) success.Sim = TRUE
      if(simulation.Attempts[i] > 100)stop("Simulation is screwing up too much")
    }
    
    #update records to include in final object
    simulated.Models[[i]] <- true.Model
    simulated.Parameters[[i]] <- pars
    simulated.Data[[i]] <- data.Vector
    if(length(train.Frame) != 0) train.Frame$y <- data.Vector
    if(length(var.Frame) != 0) var.Frame$y <- c(data.Vector, test.Vector)

    # necessary annoyance for calculating Cp later
    if("Cp" %in% names(msc.List))
      S2 <- (summary(lm(y ~ ., data = train.Frame))$sigma)^2

    # Fit all models in model.List to the newly generated data.
    fit <- function(x) try(fit.Model(fmla = x, c(data.Vector, test.Vector), var.Frame), silent=TRUE)
    fit2 <- function(x) try(fit.Model(fmla = x, data.Vector, train.Frame), silent=TRUE)
    fitted.Model.List <- lapply(model.List, fit)
    if(test.Size > 0){
      train.Model.List <- lapply(model.List, fit2)
      for(m in 1:length(model.List)) 
      {
        fitted.Model.List[[m]] = list(full = fitted.Model.List[[m]], train = train.Model.List[[m]], 
        test.Vector = test.Vector, test.Frame = test.Frame)
        if(exists("S2")) fitted.Model.List[[m]]$S2 = S2
        class(fitted.Model.List[[m]]) <- "fmo"
      }
    }
    else for(m in 1:length(model.List)) 
      {
        fitted.Model.List[[m]] = list(full = fitted.Model.List[[m]], train = NULL, test.Vector = NULL,
        test.Frame = NULL)
        if(exists("S2")) fitted.Model.List[[m]]$S2 = S2
        class(fitted.Model.List[[m]]) <- "fmo"
      }

    err.Fits <- sapply(fitted.Model.List, function(x) "try-error" %in% c(class(x$full), class(x$train)))
 
    # Sometimes fitting the true model to the simulated data gives an error... this is bad.  Just throw this iteration out.
    if(err.Fits[model.Index]==TRUE){
      if(!barebones){
        ranks.Mat[,i] <- rep(0, dim(ranks.Mat)[1])
        next}
      else {errs[i] <- 1; next}}
  
    # update the model list to include only models where the fit actually worked.
    # also change the index of the true model to reflect this updated list.
    fitted.Model.List <- fitted.Model.List[!err.Fits]
    model.Index <- which.max(sapply(model.List[!err.Fits],function(x)identical(x,true.Model)))
   
    # Whew, what a line!  Just apply all the msc functions in msc.List
    # to each fitted model, and then either take ranks of the resulting vectors
    # or standardize them, depending on whether use.Ranks is true.
    ifelse(use.Ranks, crit.Vals <- apply(t(sapply(fitted.Model.List, function(x) sapply(msc.List,
      function(y)do.call(y,list(x))))), 2, rank), crit.Vals <- apply(t(sapply(fitted.Model.List, 
      function(x) sapply(msc.List, function(y)do.call(y,list(x))))), 2, scale))
         
    # update ranks.Mat if !barebones, else update the summary stats dynamically and do not retain ranks
    # jittering breaks ties, which can make the results a bit discontinuous
    temp <- tcrossprod(crits.Mat, crit.Vals)
    temp2 <- numeric(dim(temp)[1])
    for(j in 1:nrow(temp)) {
      if(!barebones) ranks.Mat[j,i] <- sum(temp[j,] <= temp[j,model.Index])
      temp2[j] <- sum(temp[j,] <= temp[j,model.Index])
    }
    Avg.Rank <- ((i-1)/i)*Avg.Rank + (1/i)*temp2
    for(j in 1:length(thresholds))
      Probs[,j] = ((i-1)/i) * Probs[,j] + (1/i) * (temp2 > thresholds[j])
    rm(temp2)

    rm(data.Vector, true.Model, pars, fitted.Model.List, model.Index, err.Fits)
  }

  # Min rank of zero means something screwed up that iteration, so we should toss it out
  colnames(Probs) <- paste("P(Rank > ", thresholds, ")", sep='')
  colnames(crits.Mat) <- names(msc.List)
  if(!barebones) errs <- apply(ranks.Mat, 2, min)==0

  # annoying here... if sumstats is only length one, then the result of the "apply" is a vector, rather
  # than a matrix, so transposing it gives the wrong dimensions.  Hence we need these annoying if conditions.
  if(!barebones){
    if(length(sumstats)==0) Sum.Stats <- data.frame(cbind(crits.Mat, Avg.Rank, Probs))
    if(length(sumstats)==1) Sum.Stats <- data.frame(cbind(crits.Mat, apply(ranks.Mat[,!errs],1,
    function(x)sapply(sumstats,function(y)do.call(y,list(x)))), Avg.Rank, Probs))
    if(length(sumstats)>1) Sum.Stats <- data.frame(cbind(crits.Mat,  t(apply(ranks.Mat[,!errs],1,
    function(x)sapply(sumstats,function(y)do.call(y,list(x))))), Avg.Rank, Probs))
    colnames(Sum.Stats) <- c(colnames(crits.Mat), names(sumstats), "Average Rank", colnames(Probs))}
    else {Sum.Stats <- data.frame(cbind(crits.Mat, Avg.Rank, Probs))
      colnames(Sum.Stats) <- c(colnames(crits.Mat), "Average Rank", colnames(Probs))}
  if(!barebones){
    summary.Functions = as.list(rep("Internal Function", times = 1 + length(thresholds)))
    names(summary.Functions) <- c("Average Rank", colnames(Probs))
    summary.Functions <- c(sumstats, summary.Functions)
  }
  if("var.Frame" %in% search()) detach(var.Frame)
  var.Frame$y = NULL



  if(barebones){
    temp <- list(call = match.call(), Sum.Stats = Sum.Stats,
      data.Size = data.Size, num.Iter=num.Iter, msc.List=msc.List, summary.Functions = sumstats,
      var.Frame = var.Frame, model.List = model.List,
      num.Errors = sum(errs), error.Iterations = which(errs==TRUE),
      weight.Vector = weight.Vector, make.Data=deparse(substitute(make.Data)),
      stepSize=stepSize,make.Params = deparse(substitute(make.Params)),
      fit.Model = deparse(substitute(fit.Model)), simulated.Models = simulated.Models,
      simulation.Attempts = simulation.Attempts, time.Taken = proc.time() - ptm,
      thresholds = thresholds, par.Sigma = par.Sigma, data.Sigma = data.Sigma,
      test.Size = test.Size)
    class(temp) <- c( "barebones","msc")
    return(temp)
  }
  else{
    temp <- list(call = match.call(), Sum.Stats = Sum.Stats, ranks.Mat = ranks.Mat,
      data.Size = data.Size, num.Iter=num.Iter, msc.List=msc.List, summary.Functions = summary.Functions,
      var.Frame = var.Frame, model.List = model.List,
      num.Errors = sum(errs), error.Iterations = which(errs==TRUE),
      weight.Vector = weight.Vector, make.Data=deparse(substitute(make.Data)),
      stepSize=stepSize,make.Params = deparse(substitute(make.Params)),
      fit.Model = deparse(substitute(fit.Model)), simulated.Models = simulated.Models,
      simulated.Parameters = simulated.Parameters, simulated.Data = simulated.Data,
      simulation.Attempts = simulation.Attempts, time.Taken = proc.time() - ptm,
      thresholds = thresholds, par.Sigma = par.Sigma, data.Sigma = data.Sigma,
      test.Size = test.Size)

    class(temp) <- "msc"
    rm(simulated.Data, simulated.Models, simulated.Parameters)
    if(!barebones) rm(ranks.Mat)
    gc()
    return(temp)

  }
}

