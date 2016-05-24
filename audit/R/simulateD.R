simulateD <- function(ysmp,xsmp,yunsmp ,n,pgt,pwa,R)
{
  if (length(n) != 1) 
	stop("sample size n is wrong")
  if (round(n) != n || n < 1)
    stop("sample size n is not a positive integer")
  if (length(ysmp) != n)  stop("length(ysmp) != n")
  if (length(xsmp) != n)  stop("length(xsmp) != n")
  if (length(pgt) != length(pwa))
    stop("pgt and pwa not the same size")
  if (min(pwa) <= 0) stop("not all of pwa is positive")
  if (round(R) != R || R < 1)
    stop("R is not a positive integer")

	N <- n + length(yunsmp)
	obserror <-sum(ysmp-xsmp)
	obstaint <- (ysmp-xsmp)/ysmp
	obstaintnon0 <- obstaint[obstaint != 0]
	k1 <- length(obstaintnon0)
	k2 <- length(pgt)
	k <- k1+k2
	simtaints <- c(pgt,obstaintnon0,0)
	simweights <- c(pwa,rep(1,k1),n-k1)
	simD <- rep(0,R)
	
	for(i in 1:R){
		dum <- rgamma(k+1,simweights)
		dirchletparmeters <- dum/sum(dum)
		smpsize <- ceiling(dirchletparmeters[1:k]*(N-n))
		smpsize0 <- (N-n) - sum(smpsize)
		simpoptaints <- rep(simtaints,c(smpsize,smpsize0))
		smp <- sample(1:(N-n),N-n)
		simD[i] <- sum(simpoptaints*yunsmp[smp]) + obserror
	}
	return(simD)
}






