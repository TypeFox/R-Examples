#' @importFrom ape cophenetic.phylo
#' @importFrom quantreg rq
#' @importFrom vegan mantel
#' @importFrom stats lm as.dist
#' @importFrom ape cophenetic.phylo
#' @export
#' @rdname eco.xxx.regression
eco.phy.regression <- function(data,
  randomisation=c("taxa.labels", "richness", "frequency", "sample.pool", "phylogeny.pool", "independentswap", "trialswap"), permute=0, method=c("quantile", "lm", "mantel"), indep.swap=1000, abundance=TRUE, ...){
    #Assertions and argument handling
    if(! inherits(data, "comparative.comm"))  stop("'data' must be a comparative community ecology object")
    randomisation <- match.arg(randomisation)
    method <- match.arg(method)
    if(permute < 0) stop("Can't have negative null permutations!")
    
    #Setup matrices
    if(abundance==FALSE)
        data$comm[data$comm>1] <- 1
    eco.matrix <- as.dist(1 - as.matrix(comm.dist(data$comm)))
    phy.matrix <- as.dist(cophenetic.phylo(data$phy))
	
    #Observed eco.phy.regression
    observed <- .eco.phy.regression(eco.matrix, phy.matrix, method, ...)
    
    #Randomisations
    randomisations <- vector(mode="list", length=permute)
    #This won't execute if permute is 0...
    for(i in seq(from=1, length.out=permute)){
        curr.rnd <- .eco.null(data$comm, randomisation, swap.iter=indep.swap)
        rnd.mat <- as.dist(1 - as.matrix(comm.dist(curr.rnd)))
        if(any(is.na(rnd.mat))){
	    warning("NAs in permuted community matrix; skipping this iteration")
	    next()
        }
        randomisations[[i]] <- .eco.phy.regression(rnd.mat, phy.matrix, method, ...)
	}
    
    #Prepare output (...and return)
    output <- .prepare.regression.output(observed, randomisations, method, permute, "eco.phy.regression")
    output$data <- data
    return(output)
}

#Perform one set of EcoPhy regressions
.eco.phy.regression <- function(eco.mat, phy.mat, method=c("quantile", "lm", "mantel"), ...){
  method <- match.arg(method)
  if(method == 'lm')
      model <- lm(as.numeric(eco.mat) ~ as.numeric(phy.mat), ...)
  
  if(method == "quantile")
    model <- rq(as.numeric(eco.mat) ~ as.numeric(phy.mat), ...)
  
  if(method == "mantel")
    model <- mantel(eco.mat, phy.mat, ...)
	
  return(model)
}
