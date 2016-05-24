

# append res to to
# ignore first element of res since it should be the same as the last element in res (taking the last elemen as initial phi for restart)
appendCUBresults <- function(res, to)
{
  res.list.length <- length(res$b.Mat)
  init.list.length <- length(to$b.Mat)
  new.length <- res.list.length + init.list.length
  if(init.list.length == 0)
  {
    to <- list()
    to$b.Init <- res$b.Init
    to$b.RInit <- res$b.RInit
    to$p.Init <- res$p.Init
    to$phi.Init <- res$phi.Init
  }
  
  ind <- ifelse(init.list.length == 0, 1, 2)
#  s <- system.time({
  to$logL.Mat <- c(to$logL.Mat, res$logL.Mat[ind:res.list.length])
  to$b.Mat <- c(to$b.Mat, res$b.Mat[ind:res.list.length])
  to$p.Mat <- c(to$p.Mat, res$p.Mat[ind:res.list.length])
  if("phi.Mat" %in% names(res)){
    to$phi.Mat <- c(to$phi.Mat, res$phi.Mat[ind:res.list.length])
  }
  if("phi.pred.Mat" %in% names(res)){
    to$phi.pred.Mat <- c(to$phi.pred.Mat, res$phi.pred.Mat[ind:res.list.length])
  }
#  })

  return(to)
}

isConverged <- function(chains, nsamples, eps=0.1, thin=10, frac1=0.1, frac2=0.5, teston=c("phi", "sphi"), test=c("gelman", "geweke"))
{
  result <- FALSE
  dataobj <- 0
  start <- 0
  if(test[1] == "gelman")
  {
    list <- mcmc.list()
    for(i in 1:length(chains))
    {  
      if(teston[1] == "sphi")
      {
        p.mat <- do.call("rbind", chains[[i]]$p.Mat)
        index <- 3 # index of sphi for with x obs
        if(dim(p.mat)[2] == 2){index <- 2} # index of sphi for without x obs
        dataobj <- log10(p.mat[, index])
      }else if(teston[1] == "phi"){ # test on phis
        if("phi.Mat" %in% names(chains[[i]])) # with x obs 
        {
          phi.mat <- do.call("rbind", chains[[i]]$phi.Mat)
        }
        if("phi.pred.Mat" %in% names(chains[[i]])) # without x obs
        {
          phi.mat <- do.call("rbind", chains[[i]]$phi.pred.Mat)
        }
        dataobj <- log10(phi.mat)
      }else{
        stop("convergence test can not be perfomed on choosen data\n")
      }
      start <- max(0, length(dataobj) - nsamples)
      mcmcobj <- mcmc(data=dataobj, start=start, thin=thin)
      list[[i]] <- mcmcobj
    }    
    diag <- gelman.diag(list, autoburnin=FALSE)
  }else if(test[1] == "geweke"){
    if(teston[1] == "sphi")
    {
      p.mat <- do.call("rbind", chains$p.Mat)
      index <- 3 # index of sphi for with x obs
      if(dim(p.mat)[2] == 2){index <- 2} # index of sphi for without x obs
      dataobj <- log10(p.mat[, index])
    }else if(teston[1] == "phi"){ # test on phis
      if("phi.Mat" %in% names(chains)) # with x obs 
      {
        phi.mat <- do.call("rbind", chains$phi.Mat)
      }
      if("phi.pred.Mat" %in% names(chains)) # without x obs
      {
        phi.mat <- do.call("rbind", chains$phi.pred.Mat)
      }
      dataobj <- log10(phi.mat)
    }else{
      stop("convergence test can not be perfomed on choosen data\n")
    }
    start <- max(0, length(dataobj) - nsamples)
    mcmcobj <- mcmc(data=dataobj, start=start, thin=thin)
    diag <- geweke.diag(mcmcobj, frac1=frac1, frac2=frac2)
  }else{
    stop("choosen convergence test can not be found\n")
  }
  
  
  if(test[1] == "gelman")
  {
    if(teston[1] == "sphi") # scalar test on s phi
    {
      result <- abs(diag[[1]][1] - 1) < eps
      ret <- list(isConverged=result, gelman=diag[[1]][1])
    }else # else is enough here. The correctness of the method was determined above and leaves only two options
    { # multivariate test on all phi values
      result <- abs(diag$mpsrf - 1) < eps
      ret <- list(isConverged=result, gelman=diag$mpsrf)
    }
  }else { # geweke
    if(teston[1] == "sphi") # scalar test on s phi
    {
      result <- abs(diag$z) < eps
      ret <- list(isConverged=result, gelman=diag$z)
    }else{ # else is enough here. The correctness of the method was determined above and leaves only two options
      # univariate test on all phi values
      result <- sum(diag$z < eps) == length(diag$z)
      ret <- list(isConverged=result, gelman=mean(diag$z))
    }    
  }
  return(ret)
}

cubsinglechain <- function(cubmethod, frac1=0.1, frac2=0.5, reset.qr, seed=NULL, teston=c("phi", "sphi"), monitor=NULL, 
                           min=0, max=160000, conv.thin=10, eps=1, ...)
{
  ########################
  ## checking arguments ##
  ########################
  if(min > max){
    warning("min samples > max samples. setting max = min\n")
    max <- min
  }
  if(!cubmethod %in% c("cubfits", "cubappr", "cubpred")){
    stop(paste("Unkown method: ", cubmethod, "!\n"))
  }
  
  ## arguments for cub methods
  input_list <- as.list(list(...))
  
  if("p.Init" %in% names(input_list)){
    p.init <- input_list$p.Init
    input_list$p.Init <- NULL
  }else{
    p.init <- NULL
  }
  if("iterThin" %in% names(input_list)){
    # do nothing, everything is fine
  }else{
    input_list$iterThin <- 1
  }
  if("b.Init" %in% names(input_list)){
    b.init <- input_list$b.Init
    input_list$b.Init <- NULL
  }else{
    b.init <- NULL
  }  
  if("b.RInit" %in% names(input_list)){
    b.rinit <- input_list$b.RInit
    input_list$b.RInit <- NULL
  }else{
    b.rinit <- NULL
  }
  
  if("phi.Init" %in% names(input_list)){
    init.phi <- input_list$phi.Init
    input_list$phi.Init <- NULL
  }
  if("phi.pred.Init" %in% names(input_list)){
    init.pred.phi <- input_list$phi.pred.Init
    input_list$phi.pred.Init <- NULL
  }
  if(".CF.CT" %in% names(input_list)){
    .CF.CT <- input_list$.CF.CT
    input_list$.CF.CT <- NULL
  }else{
    .CF.CT <- eval(parse(text = "cubfits::.CF.CT"))
  } 
  if(".CF.CONF" %in% names(input_list)){
    .CF.CONF <- input_list$.CF.CONF
    input_list$.CF.CONF <- NULL
  }else{
    .CF.CONF <- eval(parse(text = "cubfits::.CF.CONF"))
  }   
  results <- list()
  if(is.null(seed)){
    seed <- round(runif(1, 1, 100000))
  }
  if("b.DrawScale" %in% names(input_list) ){
    b.DrawScale <- input_list$b.DrawScale
  }else{
    b.DrawScale <- .CF.CONF$b.DrawScale
  }
  if("p.DrawScale" %in% names(input_list) ){
    p.DrawScale <- input_list$p.DrawScale
  }else{
    p.DrawScale <- .CF.CONF$p.DrawScale
  }
  if("phi.pred.DrawScale" %in% names(input_list) ){
    phi.pred.DrawScale <- input_list$phi.pred.DrawScale
  }else{
    phi.pred.DrawScale <- .CF.CONF$phi.pred.DrawScale
  }
  if("phi.DrawScale" %in% names(input_list) ){
    phi.DrawScale <- input_list$phi.DrawScale
  }else{
    phi.DrawScale <- .CF.CONF$phi.DrawScale
  }   
  #############################################################
  ## running chain and checking for convergence ##
  #############################################################
  j <- 1
  gel.res <- 0
  sample.res <- 0
  converged <- FALSE
  while(!converged)
  { 
    
    
    .GlobalEnv$.CF.CT <- .CF.CT
    .GlobalEnv$.CF.CONF <- .CF.CONF
    if(cubmethod == "cubfits"){
      res <- do.call(cubfits, c(input_list, list(phi.Init = init.phi), list(p.Init = p.init), list(b.RInit = b.rinit), list(b.Init = b.init),
                                list(b.DrawScale = b.DrawScale), list(p.DrawScale = p.DrawScale), list(phi.DrawScale = phi.DrawScale)))
    }else if(cubmethod == "cubappr"){
      res <- do.call(cubappr, c(input_list, list(phi.pred.Init = init.pred.phi), list(p.Init = p.init), list(b.RInit = b.rinit), list(b.Init = b.init),
                                list(b.DrawScale = b.DrawScale), list(p.DrawScale = p.DrawScale), list(phi.pred.DrawScale = phi.pred.DrawScale)))
    }else if(cubmethod == "cubpred"){
      res <- do.call(cubpred, c(input_list, list(phi.Init = init.phi), list(phi.pred.Init = init.pred.phi), list(p.Init = p.init), list(b.RInit = b.rinit), list(b.Init = b.init),
                                list(b.DrawScale = b.DrawScale), list(p.DrawScale = p.DrawScale), list(phi.pred.DrawScale = phi.pred.DrawScale), list(phi.DrawScale = phi.DrawScale)))
    }
    b.DrawScale <- .cubfitsEnv$DrawScale$b[[length(.cubfitsEnv$DrawScale$b)]]
    if(length(.cubfitsEnv$DrawScale$p) > 0){
      p.DrawScale <- .cubfitsEnv$DrawScale$p[[length(.cubfitsEnv$DrawScale$p)]]
    }
    ## append chains and get new initial values for restart
    if(cubmethod == "cubfits" | cubmethod == "cubpred")
    {
      init.phi <- normalizeDataSet(res$phi.Mat[[length(res$phi.Mat)]])
      if(length(.cubfitsEnv$DrawScale$p) > 0){
        phi.DrawScale <- .cubfitsEnv$DrawScale$phi[[length(.cubfitsEnv$DrawScale$phi)]]
      }
    }
    if(cubmethod == "cubappr" | cubmethod == "cubpred")
    {
      init.pred.phi <- normalizeDataSet(res$phi.pred.Mat[[length(res$phi.pred.Mat)]])
      if(length(.cubfitsEnv$DrawScale$p) > 0){
        phi.pred.DrawScale <- .cubfitsEnv$DrawScale$phi.pred[[length(.cubfitsEnv$DrawScale$phi.pred)]]
      }
    }
    p.init <- res$p.Mat[[length(res$p.Mat)]]
    results <- appendCUBresults(res, results)
    
    if(length(results$p.Mat) < reset.qr) # reset the "cov" only in the begining
    {
      b.rinit <- NULL
    }else{ # use the same matrix every time after some "burnin"
      b.rinit <- res$b.RInit
    }
    currSamples <- length(results$p.Mat)
    
    ## Do convergence test
#    if(currSamples > nsamples){ #if there are not enough iterations, just keep goint until we have enough for a convergence test
      testSampleSize <- round(currSamples/2) #min(nsamples + round(growthfactor*currSamples), currSamples)
      gelman <- isConverged(results, nsamples=testSampleSize, frac1=frac1, frac2=frac2, eps=eps, thin=conv.thin, teston=teston, test="geweke")
      gel.res[j] <- gelman$gelman
      sample.res[j] <- currSamples
      cat(paste("Geweke Z score after samples: ", sample.res[j], "\t" ,gel.res[j] , "\t test was performed on ", testSampleSize/conv.thin," samples\n", sep=""))
      converged <- gelman$isConverged
      j <- j + 1
#    }
    
    ## call monitor function if it is desired
    if(!is.null(monitor))
    {
      monitor(res)
    } 
    #check if we have at least min iterations
    if(currSamples < min){converged <- FALSE}
    #check if max iteration limit is reached
    if(currSamples > max){converged <- TRUE}
  }
  ## return full length chains
  return(list(chains=results, convergence=cbind(sample.res, gel.res))) 
  
}

cubmultichain <- function(cubmethod, reset.qr, seeds=NULL, teston=c("phi", "sphi"), swap=0, swapAt=0.05, monitor=NULL, 
                          min=0, max=160000, nchains=2, conv.thin=10, eps=0.1, ncores=2, ...)
{
  #require(parallel)
  #require(doParallel)
  cl <- makeCluster(ncores)
  #registerDoParallel(cl)
  #registerDoSNOW(cl)
  
  ########################
  ## checking arguments ##
  ########################
  if(min > max){
    warning("min samples > max samples. setting max = min\n")
    max <- min
  }
  if(nchains < 2){
    stop("number of chains not sufficient. set number of chains > 2 or use cubsinglechain\n ")
  }
  if(!cubmethod %in% c("cubfits", "cubappr", "cubpred")){
    stop(paste("Unkown method: ", cubmethod, "!\n"))
  }
  
  ## arguments for cub methods
  input_list <- as.list(list(...))
  
  if("p.Init" %in% names(input_list)){
    p.init <- input_list$p.Init
    if(length(p.init) < nchains){
      stop("p.Init does not contain initial p vectors for every chain!\n")
    }
    input_list$p.Init <- NULL
  }else{
    p.init <- list(NULL)
    length(p.init) <- nchains  
  }
  if("iterThin" %in% names(input_list)){
    # do nothing, everything is fine
  }else{
    input_list$iterThin <- 1
  }  
  if("b.RInit" %in% names(input_list)){
    b.rinit <- input_list$b.RInit
    input_list$b.RInit <- NULL
  }else{
    b.rinit <- list(NULL)
    length(b.rinit) <- nchains  
  }
  if("b.Init" %in% names(input_list)){
    b.init <- input_list$b.Init
    input_list$b.Init <- NULL
  }else{
    b.init <- list(NULL)
    length(b.init) <- nchains  
  }
  
  init.phi <- list()
  init.pred.phi <- list()
  if("phi.Init" %in% names(input_list)){
    init.phi <- input_list$phi.Init
    input_list$phi.Init <- NULL
  }
  if("phi.pred.Init" %in% names(input_list)){
    init.pred.phi <- input_list$phi.pred.Init
    input_list$phi.pred.Init <- NULL
  }
  if(".CF.CT" %in% names(input_list)){
    .CF.CT <- input_list$.CF.CT
    input_list$.CF.CT <- NULL
  }else{
    .CF.CT <- eval(parse(text = "cubfits::.CF.CT"))
  } 
  if(".CF.CONF" %in% names(input_list)){
    .CF.CONF <- input_list$.CF.CONF
    input_list$.CF.CONF <- NULL
  }else{
    .CF.CONF <- eval(parse(text = "cubfits::.CF.CONF"))
  } 
  if("b.DrawScale" %in% names(input_list) ){
    b.DrawScale <- input_list$b.DrawScale
  }else{
    b.DrawScale <- list(.CF.CONF$b.DrawScale)
    length(b.DrawScale) <- nchains  
  }
  if("p.DrawScale" %in% names(input_list) ){
    p.DrawScale <- input_list$p.DrawScale
  }else{
    p.DrawScale <- list(.CF.CONF$p.DrawScale)
    length(p.DrawScale) <- nchains  
  }
  if("phi.pred.DrawScale" %in% names(input_list) ){
    phi.pred.DrawScale <- input_list$phi.pred.DrawScale
  }else{
    phi.pred.DrawScale <- list(.CF.CONF$phi.pred.DrawScale)
    length(phi.pred.DrawScale) <- nchains  
  }
  if("phi.DrawScale" %in% names(input_list) ){
    phi.DrawScale <- input_list$phi.DrawScale
  }else{
    phi.DrawScale <- list(.CF.CONF$phi.DrawScale)
    length(phi.DrawScale) <- nchains  
  } 
  
  
  
  results <- list()
  length(results) <- nchains
  if(is.null(seeds)){
    seeds <- round(runif(nchains, 1, 100000))
  }
  aa.names <- names(input_list$y)
  #############################################################
  ## running chains in parralel and checking for convergence ##
  #############################################################
  j <- 1
  gel.res <- 0
  sample.res <- 0
  converged <- FALSE
  while(!converged)
  { 
    ## run chains in parallel
    #for(i in nchains:1) # for debuging
    res <- foreach(i = 1:nchains) %dopar%
    {
      suppressMessages(library(cubfits, quietly = TRUE))
      .GlobalEnv$.CF.CT <- .CF.CT
      .GlobalEnv$.CF.CONF <- .CF.CONF
      set.seed(seeds[i])
      if(cubmethod == "cubfits"){
        do.call(cubfits, c(input_list, list(phi.Init = init.phi[[i]]), list(p.Init = p.init[[i]]), list(b.RInit = b.rinit[[i]]), list(b.Init = b.init[[i]]),
                           list(b.DrawScale = b.DrawScale[[i]]), list(p.DrawScale = p.DrawScale[[i]]), list(phi.DrawScale = phi.DrawScale[[i]])))
      }else if(cubmethod == "cubappr"){
        res <- do.call(cubappr, c(input_list, list(phi.pred.Init = init.pred.phi[[i]]), list(p.Init = p.init[[i]]), list(b.RInit = b.rinit[[i]]), list(b.Init = b.init[[i]]),
                                  list(b.DrawScale = b.DrawScale[[i]]), list(p.DrawScale = p.DrawScale[[i]]), list(phi.pred.DrawScale = phi.pred.DrawScale[[i]])))
      }else if(cubmethod == "cubpred"){
        do.call(cubpred, c(input_list, list(phi.Init = init.phi[[i]]), list(phi.pred.Init = init.pred.phi[[i]]), list(p.Init = p.init[[i]]), list(b.RInit = b.rinit[[i]]), list(b.Init = b.init[[i]]),
                           list(b.DrawScale = b.DrawScale[[i]]), list(p.DrawScale = p.DrawScale[[i]]), list(phi.pred.DrawScale = phi.pred.DrawScale[[i]]), list(phi.DrawScale = phi.DrawScale[[i]])))
      }
    }
    ## append chains and get new initial values for restart
    for(i in 1:nchains)
    {
      b.DrawScale[[i]] <- .cubfitsEnv$DrawScale$b[[length(.cubfitsEnv$DrawScale$b)]]
      p.DrawScale[[i]] <- .cubfitsEnv$DrawScale$p[[length(.cubfitsEnv$DrawScale$p)]]
      if(cubmethod == "cubfits" | cubmethod == "cubpred")
      {
        init.phi[[i]] <- normalizeDataSet(res[[i]]$phi.Mat[[length(res[[i]]$phi.Mat)]])
        phi.DrawScale[[i]] <- .cubfitsEnv$DrawScale$phi[[length(.cubfitsEnv$DrawScale$phi)]]
      }
      if(cubmethod == "cubappr" | cubmethod == "cubpred")
      {
        init.pred.phi[[i]] <- normalizeDataSet(res[[i]]$phi.pred.Mat[[length(res[[i]]$phi.pred.Mat)]])
        phi.pred.DrawScale[[i]] <- .cubfitsEnv$DrawScale$phi.pred[[length(.cubfitsEnv$DrawScale$phi.pred)]]
      }
      p.init[[i]] <- res[[i]]$p.Mat[[length(res[[i]]$p.Mat)]]
      results[[i]] <- appendCUBresults(res[[i]], results[[i]])
      
      if(length(results[[i]]$p.Mat) < reset.qr) # reset the "cov" only in the begining
      {
        cat(paste("Reset Cov matrix, total samples by now", length(results[[i]]$p.Mat), "\n"))
        b.rinit[i] <- list(NULL)
      }else{ # use the same matrix every time after some "burnin"
        b.rinit[[i]] <- res[[i]]$b.RInit
      }
    }
    currSamples <- length(results[[1]]$p.Mat)
    ## Do convergence test
    #if(currSamples > nsamples){ #if there are not enough iterations, just keep goint until we have enough for a convergence test
      testSampleSize <- round(currSamples/2) #min(nsamples + round(growthfactor*currSamples), currSamples)
      gelman <- isConverged(results, testSampleSize, eps=eps, thin=conv.thin, teston=teston, test="gelman")
      gel.res[j] <- gelman$gelman
      sample.res[j] <- currSamples
      cat(paste("Gelman score after sample: ", sample.res[j], "\t" ,gel.res[j] , "\t test was performed on ", testSampleSize/conv.thin," samples\n", sep=""))
      converged <- gelman$isConverged
      
      ## only need to swap if a convergence test was done
        ## swap bmatrix (deltat, log(mu)) 
        ## swap them here so the latest convergence check is taken into account
        #second caluse of if question not evaluated if first is false => no crash if j == 1
      if( (j > 1) && (swap > 0.0) && (swapAt > 0.0) && (abs(gel.res[j] - gel.res[j - 1]) < swapAt) )
      {
        chain.order <- sample(1:nchains, nchains)
        for(i in chain.order)
        {
          swapchain <- sample((1:nchains)[-i], 1) # avoid swaping with itself theoretically still possible since I swap in serial
          cat(paste("swapping", swap*100, "% of b matrix of chain", i, "with b matrix of chain", swapchain, "\n"))
          b.swaps <- getBInitFromSwapBMat(res[[i]]$b.Mat, res[[swapchain]]$b.Mat, swap)
          b.init[[i]] <- convert.bVec.to.b(b.swaps$b.Init1, aa.names)
          b.init[[swapchain]] <- convert.bVec.to.b(b.swaps$b.Init2, aa.names)
        }
      }else{
        b.init <- list(NULL)
        length(b.init) <- nchains  
      }
      j <- j + 1
    #}
    
    
    if(!is.null(monitor))
    {
      for(i in 1:nchains)
      {
        monitor(res, i)
      }
    } 
    
    #check if we have at least min iterations
    if(currSamples < min){converged <- FALSE}
    #check if max iteration limit is reached
    if(currSamples > max){converged <- TRUE}
  }
  ## return full length chains
  return(list(chains=results, convergence=cbind(sample.res, gel.res)))
}


getBInitFromSwapBMat <- function(b.mat1, b.mat2, swap=0.0)
{
  lastIndex <- length(b.mat1)
  numPar <- length(b.mat1[[lastIndex]])
  swapCount <- round(numPar*swap)
  
  toSwap <- sample(1:numPar, swapCount)
  temp <- b.mat1[[lastIndex]][toSwap]
  b.mat1[[lastIndex]][toSwap] <- b.mat2[[lastIndex]][toSwap]
  b.mat2[[lastIndex]][toSwap] <- temp

  return( list(b.Init1=b.mat1[[lastIndex]], b.Init2=b.mat2[[lastIndex]]) )
}







