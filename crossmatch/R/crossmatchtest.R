## crossmatch test
## z is binary vector indicating the group
## D is a distance matrix of covariates

crossmatchtest <- function(z,D)
  {
  if ( !isSymmetric(D) )
    {
    stop("Invalid distance matrix: your distance matrix is not symmetric")
    return(NA)
    }
  if ( sum(D < 0 ) > 0 )
    {
    stop("Invalid distance matrix: your distance matrix includes negative values")
    return(NA)
    }
  plainmatrix <- 100000*D/max(as.vector(D))
  diag(plainmatrix) <- 0
  mdm <- distancematrix(plainmatrix)
  nzero <- sum(mdm==0) - length(diag(mdm))
  if (nzero/((dim(mdm)[1])^2-dim(mdm)[1])>.95)
    {
    warning("Your distance matrix has some very large relative distances such that more than 95 percent of distances were rounded to zero")
    }
  res <- nonbimatch(mdm)
  mt <- pmin(as.numeric(res$matches$Group1.Row),as.numeric(res$matches$Group2.Row))
  if ( length(z) < length(mt) ) ##if the number of observations is odd remove observation that paired with ghost
    {
    mt[ mt==mt[length(mt)] ] <- 0
    mt <- mt[1:(length(mt)-1)]
    }
  z0 <- z[mt>0]
  mt0 <- factor(mt[mt>0])
  tab <- table(factor(z0),mt0)
  a1 <- sum(tab[1,]==1)
  bigN <- length(z0)
  n <- sum(z0)
  if (bigN<340) ##if the number of observations is below 340 compute the exact null distribution
    {
    dist <- crossmatchdist(bigN,n)
    pval <- dist[5,dist[2,]==a1]
    }
  else
    {
    pval<-NA
    }
  m <- bigN-n
  Ea1 <- (n*m/(bigN-1))
  Va1 <- 2*n*(n-1)*m*(m-1)/((bigN-3)*(bigN-1)*(bigN-1))
  dev <- (a1-Ea1)/sqrt(Va1)
  approx <- pnorm(dev)
  list (a1=a1,Ea1=Ea1,Va1=Va1,dev=dev,pval=pval,approxpval=approx)
  }

