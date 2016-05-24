"howmany" <-
function(pvalues,alpha=0.05,cutoff=0.05/length(pvalues),m=length(pvalues))
  {

    if(m<length(pvalues)) stop("total number of tests must cannot be smaller than number of supplied p-values")

    m0 <- length(pvalues)
    
    ord <- order(pvalues)
    pvalues <- pvalues[ord]

    howmany <- list()
    howmany$order <- ord
    howmany$pvalues <- pvalues
    howmany$alpha <- alpha
    howmany$m <- m
  
    ##cutoff
    pvalues[pvalues<cutoff] <- cutoff

    ##calculate bounding function
    boundingfunction <- get.boundingfunction.independent(m,alpha,pvalues)

    ##compute the lower bound for the number of correct rejections
    lowerbound <- numeric(m0)
    cummax <- 0
    for (p in 1:m0){

      cummax <- max(cummax,floor((p-floor(boundingfunction[p]))/max(0.2,(1-pvalues[p]))))
      lowerbound[p] <- cummax
    }

    ##give back the result
    howmany$boundingfunction <- boundingfunction
    howmany$lowerbound <- lowerbound
    class(howmany) <- "howmany"
    
    return(howmany)
  }

