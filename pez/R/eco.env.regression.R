#' @importFrom quantreg rq
#' @importFrom vegan mantel
#' @importFrom stats as.dist lm
#' @export
#' @rdname eco.xxx.regression
#' @name eco.xxx.regression
#' @importFrom stats as.dist
eco.env.regression <- function(data, randomisation=c("taxa.labels", "richness", "frequency",
                                   "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
                               permute=0, method=c("quantile", "lm", "mantel"), altogether=TRUE, indep.swap=1000, abundance=TRUE, ...){
  #Assertions and argument handling
  if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
  randomisation <- match.arg(randomisation)
  method <- match.arg(method)
  if(permute < 0) stop("Can't have negative null permutations!")
  
  #Setup matrices
  if(abundance==FALSE)
      data$comm[data$comm>1] <- 1
  eco.matrix <- as.dist(1-as.matrix(comm.dist(data$comm)))
  if(!is.null(data$env)) traits.matrix <- pianka.dist(data, altogether) else stop("'data' must contain environmental data for an environmental regression!")
  
  #Observed eco.env.regression
  if(altogether){
    observed <- .eco.env.regression(eco.matrix, traits.matrix, NULL, method, ...)} else {
      #Do separately for all traits
      observed <- vector("list", ncol(data$env))
      for(i in seq(ncol(data$env)))
        observed[[i]] <- .eco.env.regression(eco.matrix, traits.matrix, i, method, ...)
    }
  
  #Randomisations
  if(altogether){
    #Using mean of traits
    randomisations <- vector(mode="list", length=permute)
    #This won't execute if permute is 0...
    for(i in seq(from=1, length.out=permute)){
      curr.rnd <- .eco.null(data$comm, randomisation, swap.iter=indep.swap)
      rnd.mat <- as.dist(1 - as.matrix(comm.dist(curr.rnd)))
      if(any(is.na(rnd.mat))){
        warning("NAs in permuted community matrix; skipping this iteration")
        next()
      }
      randomisations[[i]] <- .eco.env.regression(rnd.mat, traits.matrix, NULL, method, ...)
    }
  } else {
    #Separately for each trait
    # - preallocate
    randomisations <- vector(mode="list", length=ncol(data$env))
    for(i in seq_along(randomisations)) randomisations[[i]] <- vector("list", permute)
    for(i in seq(from=1, length.out=permute)){
      curr.rnd <- .eco.null(data$comm, randomisation, swap.iter=indep.swap)
      rnd.mat <- as.dist(1 - as.matrix(comm.dist(curr.rnd)))
      if(any(is.na(rnd.mat))){
        warning("NAs in permuted community matrix; skipping this iteration")
        next()
      }
      for(j in seq(ncol(data$env)))
        randomisations[[j]][[i]] <- .eco.env.regression(rnd.mat, traits.matrix, j, method, ...)
    }
  }
  
  #Prepare output (...and return)
  if(altogether){
    output <- .prepare.regression.output(observed, randomisations, method, permute, "eco.env.regression")
    output$altogether <- TRUE
  } else{
      output <- vector("list", ncol(data$env))
      for(i in seq_along(output)){
        output[[i]] <- .prepare.regression.output(observed[[i]], randomisations[[i]], method, permute, "eco.env.regression")
        output[[i]]$altogether <- FALSE
      }
      output$altogether <- FALSE
      output$type <- "eco.env.regression.list"
      class(output) <- "eco.xxx.regression.list"
    }
  output$data <- data
  output$permute<-permute;output$method<-method
  return(output)
}

#Perform one set of EcoPhy regressions
.eco.env.regression <- function(eco.mat, trait.mat, which.trait=NULL,
                                method=c("quantile", "lm", "mantel"), ...){
  method <- match.arg(method)
  #Check to see if we're doing this across all traits
  if(!is.null(which.trait)) trait.mat <- as.dist(trait.mat[,,which.trait])
  
  if(method == 'lm')
    model <- lm(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "quantile")
    model <- rq(as.numeric(eco.mat) ~ as.numeric(trait.mat), ...)
  
  if(method == "mantel")
    model <- mantel(eco.mat, trait.mat, ...)
  
  return(model)
}
