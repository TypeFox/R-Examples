.statisticsHolmsAdjustemnt <- function(p){
  if (all(nna <- !is.na(p))) 
    nna <- TRUE
  p2 <- p[nna]
  n <- length(p2)
  i <- n:1
  o <- order(p2)
  ro <- order(o)
  p[nna] <- pmin(1, cummax(i* p2[o]))[ro]
  p
}

.statisticsShafferAdjustemnt <- function(p){
  n <- length(p)
  o <- order(p)
  ro <- order(o)
  # n = k * (k-1) / 2 --> solve the quadratic equation
  k <- (1+sqrt(1+8*n))/2
  t <- .shafferS(k)
  lT <- length(t)
  extT <- c(t[lT])
  # At least we are comparing 3 methods, so at least lT will be 3 ({0,1,3})
  for(i in lT:3){
    extT <- c(extT, rep(t[i-1],t[i]-t[i-1]))
  }
  pmin(1, cummax(n:1 * p[o]))[ro]
}

.shafferS <- function(k){
  prev <- list()
  for(i in 0:k){
    tmp <- .shafferSIterative(i,prev)
    #The plus 1 is because indexes of lists (must be >0)
    prev[[i+1]] <- tmp
  }
  tmp
}
.shafferSIterative <- function(k, previous){
  if(k<=1){
    mySet <- c(0)
  }
  else{
    mySet <- c()
    for(j in 1:k){
      #The plus 1 is because indexes of lists (must be >0)
      tmp <- previous[[k-j+1]]
      tmp <- tmp + choose(j,2)
      mySet <- union(mySet,tmp)
    }
  }
  mySet <- mySet[order(mySet)]
  mySet
}

.statisticsWilcoxon <- function(x,y, alternative = "less"){
  # Performs the paired wilcoxon signed rank test. Ties are allowed.
  #
  # Args:
  #   x:          a numeric array (results for the first method).
  #   y:          a numeric array (results for the second method).
  #
  # Returns:
  #   a statisticsWilcoxon object to be used in the toolbox
  
  # PARAMETER VALIDATION:
  
  #Check that x is a numeric array
  if( !is.numeric(x) )
    stop("The parameter x must be numeric.")
  #Check that y is a numeric array
  if( !is.numeric(y) )
    stop("The parameter y must be numeric.")
  #Check that x and y have the same length
  if( length(x) != length(y) )
    stop("The length of x and y must be the same.")
  #Check that alternative is equal to "less" or "greater"
  if( alternative != "less" &&  alternative != "greater")
    stop("The alternative parameter must be \"less\" or \"greater\"")
  
  N <- length(x)
  
  d <- x-y
  
  r <- rank(abs(d),ties.method = "average")
  
  #R+
  rP <- sum( r[x<y] )
  
  #R-
  rM <- sum( r[x>y] )
  
  #In case of ties
  rEq = r[x==y]
  #If it is not odd, we discard one case
  if(length(r[x==y]) %% 2==1){
    rEq <- rEq[-1]
    N <- N - 1
  }
  rEq <- sum(rEq)
  rP <- rP+rEq/2.0
  rM <- rM+rEq/2.0
  if(alternative=="less")
    t = min(rP,rM)
  else
    t = max(rP,rM)
  
  #If there is only one data in x and y, an it is the same it would be removed and N would be equal to 0.
  #In that case, psignrank(0,0) would rise a warning message.
  #Instead of that, we undo the removal of that case (only 1 draw). Therefore, 
  #we set the pvalue equal to 1 and the statistic to 0.5
  if(N==0){
    t <- 0.5
    pvalue <- 1
  }
  else
    pvalue <- psignrank(t,N)
  
  results <- list(
    "statistic"           = t,
    "pvalue"              = pvalue
  )
  class(results) <- c("statisticsWilcoxon")
  results
}

.statisticsFriedman <- function(data, rankOrder){
  # Performs the friedman test.
  #
  # Args:
  #   data:       a data.frame containing the results for each method (rows) for each problem (columns)
  #   rankOrder:  tell us if we want to maximize ("max") or minimize ("min") the output variable.
  #
  # Returns:
  #   a .statisticsFriedman object to be used in the toolbox
  
  # PARAMETER VALIDATION:
  
  #Check that variable x is a data.frame  
  if( !is.data.frame(data) || all(dim(data)==dim(data.frame())) )
    stop("The parameter data must be a non-empty data.frame")
  #Check that rankOrder is one of the two possible values
  if(rankOrder!="max" && rankOrder!="min")
    stop("The parameter rankOrder must be \"max\" or \"min\".")
  
  #Number of instantiations
  k <- dim(data)[1]
  #Number of problems
  n <- dim(data)[2]
  
  if (rankOrder == "min" )
    rankModifier <- 1
  else
    rankModifier <- -1
  
  ranks                 <- sapply(data*rankModifier, FUN=rank)
  #If k is equal to 1, sapply simplifies it to an array. In that case, we use the t function
  #which traspose x. If it is a numeric array, it just makes what we want to get (one row matrix)
  if(!is.matrix(ranks))
    ranks <- t(ranks)
  rownames(ranks)       <- rownames(data)
  friedmanRanks         <- rowMeans(ranks)
  
  #We obtain the chi-square statistic. At least k is equal to 1, so diviion by zero is avoided.
  chisqStatistic <- 12 * n / ( k * (k+1) ) * ( sum(friedmanRanks^2) - ( k * (k+1)^2 ) / 4 )
  friedmanDegreeOfFreedom     <- k-1
  friedmanPValue              <- pchisq(chisqStatistic, friedmanDegreeOfFreedom, lower.tail=FALSE)
  results <- list(
    "statistic"           = chisqStatistic,
    "degreesOfFreedom"    = friedmanDegreeOfFreedom,
    "pvalue"              = friedmanPValue,
    "n"                   = n,
    "k"                   = k,
    "ranks"               = ranks
  )
  class(results) <- c(".statisticsFriedman")
  results
}

.statisticsImanDavenport <- function(f){
  # Performs the iman davenport test.
  #
  # Args:
  #   f:          a .statisticsFriedman object
  #
  # Returns:
  #   a .statisticsImanDavenport object to be used in the toolbox
  
  #Number of instantiations
  k <- f$k
  #Number of problems
  n <- f$n
  
  if( k<2 || n < 2)
    stop("To compute the Iman Davenport test there must be at least 3 methods and at least 3 problems")
  
  chisqStatistic              <- f$statistic
  imanDavenStatistic          <- ( (n-1) * chisqStatistic ) / ( n * (k-1) - chisqStatistic )
  imanDavenDegreeOfFreedom_1  <- k-1
  imanDavenDegreeOfFreedom_2  <- (k-1) * (n-1)
  imanDavenPValue             <- pf(imanDavenStatistic, imanDavenDegreeOfFreedom_1, imanDavenDegreeOfFreedom_2, lower.tail=FALSE)
  results <- list(
    "statistic"           = imanDavenStatistic,
    "degreesOfFreedom1"   = imanDavenDegreeOfFreedom_1,
    "degreesOfFreedom2"   = imanDavenDegreeOfFreedom_2,
    "pvalue"              = imanDavenPValue,
    "n"                   = n,
    "k"                   = k,
    "ranks"               = f$ranks
  )
  class(results) <- c(".statisticsImanDavenport")
  results
}

.statisticsGeneralControl <- function(f){
  # Performs the general control test.
  #
  # Args:
  #   f:          a testFriedman object
  #
  # Returns:
  #   a .statisticsGeneralControl object to be used in the toolbox
  
  # PARAMETER VALIDATION:
  
  #Check that variable f is a testFriedman object  
  if( !is.testFriedman(f) )
    stop("The parameter f must be a testFriedman object")
  
  friedmanRanks <- rowMeans(f$ranks)
  
  #Number of instantiations
  k <- dim(f$ranks)[1]
  #Number of problems
  n <- dim(f$ranks)[2]
  
  controlIdx   <- which.min(friedmanRanks)
  criticalDiff <- sqrt(  (k*(k+1))/(6*n)   )
  
  # Generate matrix with statisctics and p-values for post-hoc test:
  p <- c()
  z <- c()
  for (i in 1:k){
    z[i]             <- abs( ( friedmanRanks[i] - friedmanRanks[controlIdx] ) / criticalDiff )
    p[i]             <- 2 * pnorm( z[i], lower.tail=FALSE)
  }
  #The p-value in the position controlIdx is equal to pnorm( 0, lower.tail=FALSE ) --> It must be NA
  p[controlIdx]      <- NA
  
  results <- list(
    "statistic"           = z,
    "criticalDifference"  = criticalDiff,
    "controlIndex"        = controlIdx,
    "pvalue"              = p
  )
  class(results) <- c(".statisticsGeneralControl")
  results
}

.statisticsGeneralPairwise <- function(f){
  # Performs the general pairwise test.
  #
  # Args:
  #   f:          a testFriedman object
  #
  # Returns:
  #   a .statisticsGeneralPairwise object to be used in the toolbox
  
  # PARAMETER VALIDATION:
  
  #Check that variable f is a testFriedman object  
  if( !is.testFriedman(f) )
    stop("The parameter f must be a testFriedman object")
  
  friedmanRanks <- rowMeans(f$ranks)
  
  # Order data according to ranks. It is necessary because we apply the test comparisson only in one direction (greater versus worse)
  correct_data_order <- order(friedmanRanks)
  friedmanRanks      <- friedmanRanks[correct_data_order]
  algorithmNames     <- names(friedmanRanks)
  
  #Number of instantiations
  k <- dim(f$ranks)[1]
  #Number of problems
  n <- dim(f$ranks)[2]
  
  criticalDiff <- sqrt(  (k*(k+1))/(6*n)   )
  
  met1       <- c()
  met2       <- c()
  p          <- c()
  z          <- c()
  counter    <- 1
  for (i in 1:(k-1))
  {
    for (j in (i+1):k)
    {
      met1[counter]         <- algorithmNames[i]
      met2[counter]         <- algorithmNames[j]
      z[counter]            <- abs( ( friedmanRanks[j] - friedmanRanks[i] ) / criticalDiff )
      p[counter]            <- 2 * pnorm( z[counter], lower.tail=FALSE)
      counter               <- counter + 1
    }
  }
  
  names <- data.frame( method1 = met1, method2 = met2)
    
  results <- list(
    "statistic"           = z,
    "criticalDifference"  = criticalDiff,
    "names"               = names,
    "pvalue"              = p
  )
  class(results) <- c(".statisticsGeneralPairwise")
  results
}
