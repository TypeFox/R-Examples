SteelConfInt <-function(..., data = NULL, conf.level=.95,
		alternative=c("less","greater","two.sided"),
		method = c("asymptotic","exact","simulated"), Nsim = 10000){

pmaxWilcox <- function(x,ns){
#----------------------------------------------------------------------------------
# computes the CDF of the maximum of standarized Wilcoxon statistics for each 
# element in x, using the asymptotic approximation.
# Standardized Wilcoxon statistics are considered for sample sizes ns[1] and 
# each of ns[i], i=2,..., length(ns)
#----------------------------------------------------------------------------------
	if(length(ns) < 2) return("ns has length < 2\n")
	n1 <- ns[1]
	ni <- ns[-1]
	s <- length(ni)
	f1 <- sqrt(1+ni/(n1+1))
	f2 <- sqrt(ni/(n1+1))
	fx <- function(z,U,f1,f2){
 		fxx <- 1
		for(i in 1:s){
			fxx <- fxx * pnorm(U*f1[i]-z*f2[i])
		}
		fxx*dnorm(z)
	}
	N <- length(x)
	px <- numeric(N)
	for(j in 1:N){
		px[j] <- integrate(f=fx,-Inf,Inf,U=x[j],f1=f1,f2=f2)$value
	}
px
}
# end of pmaxWilcox
qmaxWilcox <- function(p,ns){
#----------------------------------------------------------
# computes the p-quantile of the CDF provided by pmaxWilcox
# not vectorized over p
#----------------------------------------------------------
	if(length(ns) < 2) return("ns has length < 2\n")
	fx <- function(z,ns,p){
 		pmaxWilcox(z,ns)-p
	}
	a <- qnorm(p)-2
	b <- qnorm(p)+2
	while(fx(a,ns,p) > 0) a <- a-1
	while(fx(b,ns,p) < 0) b <- b+1
	uniroot(fx,c(a,b),ns,p)$root
}
# end of qmaxWilcox
ProbWilcox <- function(U,n1,ni){
#-----------------------------------------------------------------------
# computes the asymptotic joint probability P(W.i <= U[i], i=1,\ldots,s)
# where W.i is the Mann-Whitney statistic for sample sizes n1 and ni[i], 
# for i=1,...,s, and the same sample is involved for the first sample
# of size n1 in all Mann-Whitney statistics.
#-----------------------------------------------------------------------
	s <- length(ni)
	if(s < 1) return("ni has length < 1\n")
        if(length(U) != s) return(paste("U does not have length",s,"\n"))
	f1 <- sqrt(1+ni/(n1+1))
	f2 <- sqrt(ni/(n1+1))
	fx <- function(z,U,f1,f2){
 		fxx <- 1
		for(i in 1:s){
			fxx <- fxx * pnorm(U[i]*f1[i]-z*f2[i])
		}
		fxx*dnorm(z)
	}
	px <- integrate(f=fx,-Inf,Inf,U=U,f1=f1,f2=f2)$value
	px
}
# end of ProbWilcox
qdiscrete <- function (x, gam) 
{
#------------------------------------------------------------
# x is assumed to be a sorted numeric vector and 0 < gam < 1
# This function finds cm and cp such that 
# cm is the largest value such that the proportion of 
# x values <= cm is <= gam. cm is set to -Inf, when the 
# highest achievable proportion <= gam is zero.
# cp is the smallest value such that the proportion of 
# x values <= cp is >= gam
#------------------------------------------------------------
    if (gam <= 0 | gam >= 1) 
        stop(paste("conf.level not in (0,1)","\n"))
    N <- length(x)
    gamN <- gam * N
    im <- floor(gamN)
    ip <- ceiling(gamN)
    cm <- -Inf
    cp <- x[ip]
    if(im > 0){ 
	if(x[im] < x[ip]){ cm <- x[im]}else{ 
		if(sum(x <= x[im]) <= gamN){cm <- x[im]}else{
		 	if(x[im] > x[1]) cm <- max(x[x<x[im]])
		}
	}
    }
    list(cm = cm, cp = cp)
}
# end of qdiscrete function
samples <- io(..., data = data)
alternative <- match.arg(alternative)
method <- match.arg(method)
out <- na.remove(samples)
na.t <- out$na.total
if( na.t > 1) cat(paste("\n",na.t," NAs were removed!\n\n"))
if( na.t == 1) cat(paste("\n",na.t," NA was removed!\n\n"))
samples <- out$x.new
k <- length(samples)
if (k < 2) stop("Must have at least two samples.")
ns <- sapply(samples, length)
nsamp <- sum(ns)
if (any(ns == 0)) stop("One or more samples have no observations.")
x <- numeric(nsamp)
istart <- 0
for (i in 1:k){
	x[istart+(1:ns[i])] <- samples[[i]]
	istart <- istart + ns[i]		
}
if(method != "asymptotic"){
	if(!is.loaded("SteelConf")) dyn.load("SteelConf.so")
}
if(alternative=="two.sided"){gam <- 1-(1-conf.level)/2 }else{gam <- conf.level}
if(length(x) != sum(ns)) return("sample sizes don't agree with length of data vector\n")
k <- length(ns)
s <- k-1
# number of treatments
rx <- rank(x)
n1 <- ns[1]
ni <- ns[-1]
mu <- ni*n1/2 
	# vector of means for Mann-Whitney statistics
tau <- sqrt(n1*ni*(n1+ni+1)/12) 
	# vector of stand. deviations of Mann-Whitney statistics
	# in the case of continuous sampled populations (no ties).  
n.ties <- nsamp - length(unique(x))
# computing total number of combination splits
ncomb <- choose(nsamp,ns[1])
np <- nsamp-ns[1]
if(k>2){
	for(i in 2:(k-1)){
		ncomb <- ncomb * choose(np,ns[i])
        	np <- np-ns[i]
	}
}
if(method == "exact"){
	useExact <- TRUE}else{
        useExact <- FALSE
}
# calculation of coverage probabilities based on asymptotics
cgam <- qmaxWilcox(gam,ns)
if(alternative != "greater"){
	ellpU <- ceiling(tau*cgam+mu+1)
	ellpUx <- ellpU
        ellpUx[ellpUx > n1*ni] <- Inf
	ellmU <- floor(tau*cgam+mu+1)
	ellmUx <- ellmU
        ellmUx[ellmUx > n1*ni] <- Inf
	ellcU <- round(tau*cgam+mu+1)
	ellcUx <- ellcU
        ellcUx[ellcUx > n1*ni] <- Inf
	probpU <- ProbWilcox((ellpUx-1-mu)/tau,n1,ni)
	probmU <- ProbWilcox((ellmUx-1-mu)/tau,n1,ni)
	probcU <- ProbWilcox((ellcUx-1-mu)/tau,n1,ni)
}
if(alternative != "less"){
	ellmL <- ceiling(n1*ni-tau*cgam-mu)
	ellmLx <- ellmL
	ellmLx[ellmLx < 1] <- -Inf
	ellpL <- floor(n1*ni-tau*cgam-mu)
	ellpLx <- ellpL
	ellpLx[ellpLx < 1] <- -Inf
	ellcL <- round(n1*ni-tau*cgam-mu)
	ellcLx <- ellcL
	ellcLx[ellcLx < 1] <- -Inf
	probpL <- ProbWilcox((n1*ni-ellpLx-mu)/tau,n1,ni)
	probmL <- ProbWilcox((n1*ni-ellmLx-mu)/tau,n1,ni)
	probcL <- ProbWilcox((n1*ni-ellcLx-mu)/tau,n1,ni)
}
# initializing asymptotic bounds
Lboundp <- rep(-Inf,s)
Uboundp <- rep(Inf,s)
Lboundc <- rep(-Inf,s)
Uboundc <- rep(Inf,s)
# initializing bounds based on exact or simulated null distribution
LboundXp <- rep(-Inf,s)
UboundXp <- rep(Inf,s)
LboundXc <- rep(-Inf,s)
UboundXc <- rep(Inf,s)

if(alternative != "greater"){
# determine ellUp and probUp for which probUp comes closest 
# to gam but >= gam among the three possibilities ellpU,
# ellcU, and ellmU, with corresponding probpU, probcU, probmU
	ellUp <- ellpU
	probUp <- probpU
        if(probcU >= gam & probcU < probUp){
		ellUp <- ellcU
                probUp <- probcU
	}
        if(probmU >= gam & probmU < probUp){
		ellUp <- ellmU
                probUp <- probmU
	}
	ellUc <- ellcU
	probUc <- probcU
# determine ellUc and probUc for which probUc comes closest 
# to gam among the three possibilities ellpU,
# ellcU, and ellmU, with corresponding probpU, probcU, probmU
	if(abs(probUc-gam)>abs(probpU-gam)){
		ellUc <- ellpU
		probUc <- probpU
	}
	if(abs(probUc-gam)>abs(probmU-gam)){
		ellUc <- ellmU
		probUc <- probmU
	}

}
if(alternative != "less"){
# determine ellLp and probLp for which probLp comes closest 
# to gam but >= gam among the three possibilities ellpL,
# ellcL, and ellmL, with corresponding probpL, probcL, probmL
		ellLp <- ellpL
		probLp <- probpL
           	if(probcL >= gam & probcL < probLp){
			ellLp <- ellcL
                 	probLp <- probcL
		}
           	if(probmL >= gam & probmL < probLp){
			ellLp <- ellmL
                 	probLp <- probmL
		}
# determine ellLc and probLc for which probLc comes closest 
# to gam among the three possibilities ellpL,
# ellcL, and ellmL, with corresponding probpL, probcL, probmL
		ellLc <- ellcL
		probLc <- probcL
		if(abs(probLc-gam)>abs(probpL-gam)){
			ellLc <- ellpL
			probLc <- probpL
		}
		if(abs(probLc-gam)>abs(probmL-gam)){
			ellLc <- ellmL
			probLc <- probmL
		}
}
if(alternative == "less"){
	achieved.confidence.p <- probUp
	achieved.confidence.c <- probUc
}
if(alternative == "greater"){
	achieved.confidence.p <- probLp
	achieved.confidence.c <- probLc
}
if(alternative == "two.sided"){
	achieved.confidence.p <- probUp+probLp-1
	achieved.confidence.c <- probUc+probLc-1
}


# end of asymptotic treatment
# calculation based on exact enumeration or simulation
if(method != "asymptotic"){
	if((Nsim < ncomb & method == "exact") |  method == "simulated"){
		method <- "simulated"
		useExact <- FALSE
	   	nrow <- Nsim * s
	}else{
		nrow <- ncomb * s
	}
	MannWhitneyStats <- numeric(nrow)

	out <- .C("SteelConf", Nsim=as.integer(Nsim), k=as.integer(k), 
		rx=as.double(rx), ns=as.integer(ns), 
            	useExact=as.integer(useExact),
		MannWhitneyStats = as.double(MannWhitneyStats))

	MannWhitneyStats <- matrix(out$MannWhitneyStats,ncol=s,byrow=T)
	maxstandMannWhitneyStats <- apply(matrix(
		(out$MannWhitneyStats-mu)/tau,ncol=s,byrow=T),1,max)
	if(method=="simulated"){
# taking advantage (only in the simulated case) of the equality of the full
# distributions of 
# -minstandMannWhitneyStats and maxstandMannWhitneyStats. 
	minstandMannWhitneyStats <- apply(matrix(
		(out$MannWhitneyStats-mu)/tau,ncol=s,byrow=T),1,min)
        maxstandMannWhitneyStats <- sort(c(maxstandMannWhitneyStats,
		-minstandMannWhitneyStats))
	}else{
		maxstandMannWhitneyStats <- sort(maxstandMannWhitneyStats)
	}
	cmp <- qdiscrete(maxstandMannWhitneyStats,gam)
	cm <- cmp$cm
	cp <- cmp$cp
	ellhat <- ceiling(tau*cp+mu)
	elltilde <- floor(tau*cm+mu)
	ellest <- round(tau*(cm+cp)/2+mu)
	selhat <- MannWhitneyStats[,1] <= ellhat[1]
	seltilde <- MannWhitneyStats[,1] <= elltilde[1]
	selest <- MannWhitneyStats[,1] <= ellest[1]
	if(s >1){
		for(i in 2:s){
			selhat <- selhat & (MannWhitneyStats[,i] <= ellhat[i])
			seltilde <- seltilde & (MannWhitneyStats[,i] <= elltilde[i])
			selest <- selest & (MannWhitneyStats[,i] <= ellest[i])
		}
	}
	p.ellhat <- mean(selhat)
	p.elltilde <- mean(seltilde)
	p.ellest <- mean(selest)
	if(alternative == "less"){
		ellhatU <- ellhat+1
		elltildeU <- elltilde+1
		ellestU <- ellest+1
	}
	if(alternative == "greater"){
		ellhatL <- n1*ni-ellhat
		elltildeL <- n1*ni-elltilde
		ellestL <- n1*ni-ellest
	}
	if(alternative == "two.sided"){
		ellhatU <- ellhat+1
		elltildeU <- elltilde+1
		ellestU <- ellest+1
		ellhatL <- n1*ni-ellhat
		elltildeL <- n1*ni-elltilde
		ellestL <- n1*ni-ellest
		p.ellhat <- 1-2*(1-p.ellhat)
		p.elltilde <- 1-2*(1-p.elltilde)
		p.ellest <- 1-2*(1-p.ellest)
# reset gam to original confidence level for later 
# refinement comparisons.
		gam <- conf.level
	}
	if(alternative != "greater"){
# determine ellUXp and probXp for which probXp comes closest 
# to gam but >= gam among the three possibilities ellhatU,
# elltildeU, and ellestU, with corresponding p.ellhat, p.elltilde, p.ellest
		ellUXp <- ellhatU
		probXp <- p.ellhat
        	if(p.elltilde >= gam & p.elltilde < probXp){
			ellUXp <- elltildeU
        	        probXp <- p.elltilde
		}
        	if(p.ellest >= gam & p.ellest < probXp){
			ellUXp <- ellestU
        	        probXp <- p.ellest
		}
		ellUXc <- ellhatU
		probXc <- p.ellhat
# determine ellUXc and probXc for which probXc comes closest 
# to gam among the three possibilities ellhatU,
# elltildeU, and ellestU, with corresponding p.ellhat, p.elltilde, p.ellest
		if(abs(probXc-gam)>abs(p.elltilde-gam)){
			ellUXc <- elltildeU
			probXc <- p.elltilde
		}
		if(abs(probXc-gam)>abs(p.ellest-gam)){
			ellUXc <- ellestU
			probXc <- p.ellest
		}
	}
	if(alternative != "less"){
# determine ellLXp and probXp for which probXp comes closest 
# to gam but >= gam among the three possibilities ellhatL,
# elltildeL, and ellestL, with corresponding p.ellhat, p.elltilde, p.ellest
		ellLXp <- ellhatL
		probXp <- p.ellhat
	        if(p.elltilde >= gam & p.elltilde < probXp){
			ellLXp <- elltildeL
	                probXp <- p.elltilde
		}
	        if(p.ellest >= gam & p.ellest < probXp){
			ellLXp <- ellestL
	                probXp <- p.ellest
		}
		ellLXc <- ellhatL
		probXc <- p.ellhat
# determine ellLXc and probXc for which probXc comes closest 
# to gam among the three possibilities ellhatL,
# elltildeL, and ellestL, with corresponding p.ellhat, p.elltilde, p.ellest
		if(abs(probXc-gam)>abs(p.elltilde-gam)){
			ellLXc <- elltildeL
			probXc <- p.elltilde
		}
		if(abs(probXc-gam)>abs(p.ellest-gam)){
			ellLXc <- ellestL
			probXc <- p.ellest
		}
	}
}
if(alternative != "greater"){
	istart <- n1
	for( i in 1:s ){
		D <- sort(outer(x[(istart+1):(istart+ni[i])],x[1:n1],"-"))
		if(ellUp[i] > n1*ni[i]){ 
			Uboundp[i] <- Inf}else{
			Uboundp[i] <- D[ellUp[i]]
		}
		if(ellUc[i] > n1*ni[i]){ 
			Uboundc[i] <- Inf}else{
			Uboundc[i] <- D[ellUc[i]]
		}
		if(method != "asymptotic"){
			if(ellUXp[i] > n1*ni[i]){ 
				UboundXp[i] <- Inf}else{
				UboundXp[i] <- D[ellUXp[i]]
			}
			if(ellUXc[i] > n1*ni[i]){ 
				UboundXc[i] <- Inf}else{
				UboundXc[i] <- D[ellUXc[i]]
			}
		}
		istart <- istart + ni[i]
	}
}
if(alternative != "less"){
	istart <- n1
	for( i in 1:s ){
		D <- sort(outer(x[(istart+1):(istart+ni[i])],x[1:n1],"-"))
		if(ellLp[i] < 1){
			Lboundp[i] <- -Inf}else{	
			Lboundp[i] <- D[ellLp[i]]
		}
		if(ellLc[i] < 1){
			Lboundc[i] <- -Inf}else{	
			Lboundc[i] <- D[ellLc[i]]
		}
		if(method != "asymptotic"){
			if(ellLXp[i] < 1){
				LboundXp[i] <- -Inf}else{	
				LboundXp[i] <- D[ellLXp[i]]
			}
			if(ellLXc[i] < 1){
				LboundXc[i] <- -Inf}else{	
				LboundXc[i] <- D[ellLXc[i]]
			}
		}
		istart <- istart + ni[i]
	}
}


Delta <- paste("Delta_",1:s,sep="")
LBound <- Lboundp
UBound <- Uboundp
LBoundc <- Lboundc
UBoundc <- Uboundc
levelA <- c(round(achieved.confidence.p,6),rep("",length(LBound)-1))
levelB <- c(round(achieved.confidence.c,6),rep("",length(LBound)-1))

outA<- data.frame(L=LBound,U=UBound,level=levelA,row.names=Delta)
outB<- data.frame(L=LBoundc,U=UBoundc,level=levelB,row.names=Delta)
if(alternative=="greater"){
	ellUc <- NA
	ellUp <- NA
	ellUXc <- NA
	ellUXp <- NA
}
if(alternative=="less"){
	ellLc <- NA
	ellLp <- NA
	ellLXc <- NA
	ellLXp <- NA
}
i.LU <- cbind(ellLp,ellUp,ellLc,ellUc)
dimnames(i.LU) <- list(NULL,c("i.L","i.U","i.Lc","i.Uc"))
i.LUX <- NA
if(method=="asymptotic"){
out <- list(conservative.bounds.asymptotic = outA,
	closest.bounds.asymptotic=outB)
}else{
	i.LUX <- cbind(ellLXp,ellUXp,ellLXc,ellUXc)
	dimnames(i.LUX) <- list(NULL,c("i.L","i.U","i.Lc","i.Uc"))
	levelAX <- c(round(probXp,6),rep("",length(LboundXp)-1))
	levelBX <- c(round(probXc,6),rep("",length(LboundXp)-1))
	outAX <- data.frame(L=LboundXp,U=UboundXp,level=levelAX,row.names=Delta)
	outBX <- data.frame(L=LboundXc,U=UboundXc,level=levelBX,row.names=Delta)

	if(method == "simulated"){
		
		out <- list(conservative.bounds.asymptotic = outA, 
			closest.bounds.asymptotic=outB,
			conservative.bounds.simulated = outAX, 
			closest.bounds.simulated=outBX)
	}
	if(method == "exact"){
		out <- list(conservative.bounds.asymptotic = outA, 
			closest.bounds.asymptotic=outB,
			conservative.bounds.exact = outAX, 
			closest.bounds.exact=outBX)
	}
}

outx <- list(test.name = "Steel.bounds",
	 	n1 = n1, ns = ni, N = nsamp,n.ties=n.ties,
		bounds=out,method=method,Nsim=Nsim,i.LU=i.LU,i.LUX=i.LUX)
class(outx) <- "kSamples"
outx		
}
