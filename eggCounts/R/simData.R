simData1s <- function(n=10, mean=500, kappa=0.5, phi=1, f=50,rounding=TRUE){
  # error checking
  if (n%%1!=0) {stop("sample size is not an integer")}
  if (kappa<0) {stop("overdispersion parameter is not positive")}
  if ((phi<0)||(phi>1)) {stop("prevalence is not between 0 and 1")}
  if (sum(f%%1!=0)>0) {stop("correction factor is not integer(s)")}
  if (!is.logical(rounding)) {stop("rounding argument is not a logical")}
  if (length(f)>1 && length(f)!=n) {stop("the length of correction factors and sample size do not match")}
   # sample true egg counts
   trueMean <- rgamma(n,shape=kappa, rate=kappa/mean)
   infected <- sample(c(FALSE,TRUE),size=n, replace=TRUE,prob=c(1-phi,phi))
   trueMean[!infected] <- 0
   trueEPG <- rpois(n,trueMean)
   # take a subsample and count eggs y, such that epg = f*y
   if (rounding){
   fec <- rpois(n,lambda=round(trueEPG/f))
   } else {
     fec <- rpois(n,lambda=trueEPG/f)}
   # fec <- rbinom(n,trueEPG,1/f)
   data <- cbind(obs=fec*f, master=fec, true=trueEPG)
   return(data)
}

simData2s <- function(n=10, preMean=500, delta=0.1, kappa=0.5,phiPre=1,phiPost=phiPre, f=50, paired=TRUE,rounding=TRUE){
  # error checking
  if (n%%1!=0) {stop("sample size is not an integer")}
  if (kappa<0) {stop("overdispersion parameter is not positive")}
  if ((phiPre<0)||(phiPre>1)) {stop("pre-treatment prevalence is not between 0 and 1")}
  if ((phiPost<0)||(phiPost>1)) {stop("post-treatment prevalence is not between 0 and 1")}
  if (sum(f%%1!=0)>0) {stop("correction factor is not integer(s)")}
  if (!is.logical(paired)) {stop("paired is not a logical")}
  if (!is.logical(rounding)) {stop("rounding argument is not a logical")}
  if (length(f)>1 && length(f)!=n) {stop("the length of correction factors and sample size do not match")}
  if(paired){
	# sample true egg counts before treatment
	truePreMean <- rgamma(n,shape=kappa, rate=kappa/preMean)
	infected <- sample(c(FALSE,TRUE),size=n, replace=TRUE,prob=c(1-phiPre,phiPre))
	truePreMean[!infected] <- 0
	truePreEPG <- rpois(n,truePreMean)
	# take a subsample and count eggs y, such that epg = f*y
	if (rounding){
	fec <- rpois(n,lambda=round(truePreEPG/f))
	} else {fec <- rpois(n,lambda=truePreEPG/f)}
	# fec <- rbinom(n,truePreEPG,1/f)	
	# now add sample egg counts after treatment
	infectedA <- infected
	infectedA[infected] <- sample(c(FALSE,TRUE),size=sum(infected), replace=TRUE,prob=c(1-phiPost/phiPre,phiPost/phiPre))

	truePostMean <- truePreMean*delta
	truePostMean[!infectedA] <- 0
	truePostEPG <- rpois(n,truePostMean)
	# take a subsample and count eggs y, such that epg = f*y
	if (rounding){
	postFEC <- rpois(n,lambda=round(truePostEPG/f))
	} else {postFEC <- rpois(n,lambda=truePostEPG/f)}
	#postFEC <- rbinom(n,truePostEPG,1/f)	
	data <- cbind(obsPre =fec*f, masterPre=fec, truePre=truePreEPG, obsPost=postFEC*f, masterPost=postFEC, truePost=truePostEPG)

  } else {
	dataPre <- simData1s(n=n, mean=preMean, kappa=kappa, phi=phiPre,rounding=rounding)
	dataPost <- simData1s(n=n, mean=preMean*delta, kappa=kappa, phi=phiPost,rounding=rounding)

	data <- cbind(dataPre,dataPost)
    colnames(data) <- paste(colnames(data),rep(c("pre","post"),each=3),sep=".")
	  
  }
   return(data)
}

