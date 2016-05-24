"howmany_dependent" <-
function(X,Y,alpha=0.05,test=wilcox.test,alternative="two.sided",n.permutation=round(20/alpha))
  {
   
    if(length(dim(X))!=2) stop(" dimension of X must have length 2\n ")
    if(nrow(X)!=length(Y)) stop(" number of rows of X must match length of Y\n ")
    ##check binary class variable Y
    if(is.factor(Y)) Y <- as.numeric(Y)
    if(length(unique(Y))!=2) stop(" Y must have binary values\n ")
    valY <- unique(Y)
    Ynew <- numeric(length(Y))
    Ynew[Y==min(valY)] <- 0
    Ynew[Y==max(valY)] <- 1
    Y <- Ynew
    
    if(ncol(X)<2) stop( " number of columns of X must exceed 1\n " )
        
    m <- ncol(X)
    n <- nrow(X)

    ##calculate original p-values
    pvalues <- numeric(m)
    for (p in 1:m){
      pvalues[p] <- test(X[Y==0,p],X[Y==1,p],alternative=alternative)$p.value
    }

    ##order p-values
    ord <- order(pvalues)
    pvalues <- pvalues[ord]

    ##compute bounding function
    boundingfunction <- get.boundingfunction.dependent(X,Y,alpha,test,alternative,at=pvalues,n.permutation=max(100,n.permutation))

    ##compute lower bound for the number of correct rejections
    lowerbound <- numeric(m)
    cummax <- 0
    for (p in 1:m){
      cummax <- max(cummax,(p-floor(boundingfunction[p])))
      lowerbound[p] <- cummax
    }

    ##give back result
    howmany <- list()
    howmany$order <- ord
    howmany$pvalues <- pvalues
    howmany$alpha <- alpha
    howmany$boundingfunction <- boundingfunction
    howmany$lowerbound <- lowerbound
    howmany$m <- length(pvalues)
    class(howmany) <- "howmany"
    return(howmany)
  }

