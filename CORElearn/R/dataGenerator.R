classDataGen <- function(noInst, t1=0.7, t2=0.9, t3=0.34, t4=0.32, p1=0.5, classNoise=0)
{
	stopifnot(0 <= t1 & t1 <= 1)
	stopifnot(0 <= t2 & t2 <= 1)
	stopifnot(0 <= t3 & t3 <= 1)
	stopifnot(0 <= t4 & t4 <= 1)
	stopifnot(0 <= p1 & p1 <= 1)
	prob <- c(t1, t1, t2, t2, t2, t2, t3, t3 * (2 - t3), t3, p1, classNoise)
	mult <- c(3, 3, 6)
	n <- length(prob) + length(mult) + 3
	r <- matrix(runif(n*noInst), nrow=noInst, ncol=n, byrow=TRUE)
	for (i in 1:length(prob)) {
		r[,i] <- r[,i] < prob[i]
	}
	for (i in 11+(1:length(mult))) {
		r[,i] <- ceiling(mult[i-11] * r[,i])
	}
	for (i in 11+(4:6)) {
		r[,i] <- qnorm(r[,i])
	}
	a <- matrix(nrow=noInst, ncol=7)
	# discrete attributes for class 1
	a[,1] <- r[,1] * r[,3]
	a[,2] <- (r[,1] * r[,4] + r[,2] * r[,5]) >= 1
	a[,3] <- 1 + r[,2] * r[,6] * r[,11+1]
	a[,4] <- r[,6+1]
	a[,5] <- r[,6+2]
	a[,6] <- 1 + r[,6+3] * r[,11+2]
	a[,7] <- pmax(1, r[,11+3] - 2)
	# continuous attribute specific for class 1
	r[,11+5] <- r[,11+5] * t4
	# case classes
	cl <- 2 - r[,10]
	ind <- which(cl == 2)
	if (length(ind) >= 1) {
		a[ind,1:6] <- a[ind,c(4:6,1:3),drop=FALSE]
		r[ind,11+(4:5)] <- r[ind,11+(5:4),drop=FALSE]
	}
	# class noise
	ind <- which(r[,11] == 1)
	if (length(ind) >= 1) {
		cl[ind] <- 3 - cl[ind]
	}
	# output
	data.frame(
		a1=factor(a[,1], levels=c(0,1)),
		a2=factor(a[,2], levels=c(0,1)),
		a3=factor(letters[a[,3]], levels=c("a","b","c","d")),
		a4=factor(a[,4], levels=c(0,1)),
		a5=factor(a[,5], levels=c(0,1)),
		a6=factor(letters[a[,6]], levels=c("a","b","c","d")),
		a7=factor(letters[a[,7]], levels=c("a","b","c","d")),
		x1=r[,11+4],
		x2=r[,11+5],
		x3=r[,11+6],
		class=factor(cl, levels=c(1,2)))
}

regDataGen <- function(noInst, t1=0.8, t2=0.5, noise=0.1)
{
	stopifnot(0 <= t1 & t1 <= 1)
	prob <- c(1/2, t1, t1, 1/2, 1/2, 1/2, 1/2)
	n <- length(prob) + 9
	r <- matrix(runif(n*noInst), nrow=noInst, ncol=n, byrow=TRUE)
	for (i in 1:length(prob)) {
		r[,i] <- r[,i] < prob[i]
	}
	for (i in 7+(1:9)) {
		r[,i] <- qnorm(r[,i])
	}
	a <- matrix(nrow=noInst, ncol=7)
	# discrete attributes for r[,1] = 0
	a[,1] <- r[,2]
	a[,2] <- 2*r[,3] + r[,4] + 1
	a[,3] <- r[,5]
	a[,4] <- 2*r[,6] + r[,7] + 1
	# auxiliary continuous variables
	x1 <- r[,7+1]
	x2 <- r[,7+2]
	x3 <- r[,7+3]
	x4 <- 1/(1 + exp( - r[,7+4]))
	x5 <- 1/(1 + exp( - r[,7+5]))
	x6 <- 1/(1 + exp( - r[,7+6]))
	# internal split
    f <- numeric(noInst)
	linear <- r[,1] == 0
	if (any(linear)) {
		f[linear] <- x4[linear] - 2*x5[linear] + 3*x6[linear]
		x1[linear] <- t2 * x1[linear]
	}
	nonlinear <- !linear
	if (any(nonlinear)) {
		f[nonlinear] <- cos(4*pi*x4[nonlinear])*(2*x5[nonlinear]-3*x6[nonlinear])
		a[nonlinear,1] <- 1 - a[nonlinear,1]
		a[nonlinear,2] <- 5 - a[nonlinear,2]
		x2[nonlinear] <- t2 * x2[nonlinear]
	}
	# output
	data.frame(
		a1=factor(a[,1], levels=c(0,1)),
		a2=factor(letters[a[,2]], levels=c("a","b","c","d")),
		a3=factor(a[,3], levels=c(0,1)),
		a4=factor(letters[a[,4]], levels=c("a","b","c","d")),
		x1=x1,
		x2=x2,
		x3=x3,
		x4=1/(1 + exp( - r[,7+4] + noise * r[,7+7])),
		x5=1/(1 + exp( - r[,7+5] + noise * r[,7+8])),
		x6=1/(1 + exp( - r[,7+6] + noise * r[,7+9])),
		response=f)
}

# generate ordinal data to be used for example with ordEval algorithm
ordDataGen<-function(noInst, classNoise=0) {
    len <- noInst
    maxValue <- 5
    # generate performance, basic, random, and excitement attributes
    P <- list()
    P[[1]]<-distGen(len,c(0.2,0.2,0.2,0.2,0.2))
    P[[2]]<-distGen(len,c(0.2,0.2,0.2,0.2,0.2))
    B<-list()
    for (i in 1:2) {
        B[[i]] <- distGen(len,c(0.2,0.2,0.2,0.2,0.2))
    }
    E <- list()
    for (i in 1:2)
        E[[i]]<-distGen(len,c(0.2,0.2,0.2,0.2,0.2))
    R<-list()
    R[[1]] <- distGen(len,c(0.2,0.2,0.2,0.2,0.2))
    rn<-rnorm(len)
    R[[2]]<-1+as.integer(((rn-min(rn))/(max(rn)-min(rn)+1e-10))*5)
    # generate response
    C=vector(mode="integer",length=len)
    C<- performanceWeak(P[[1]])+performanceStrong(P[[2]])+basicWeak(B[[1]])+basicStrong(B[[2]])+
            excitementWeak(E[[1]])+excitementStrong(E[[2]])
    Cnorm <- forceDist(C,c(0.2,0.2,0.20,0.2,0.2))
    #class noise
    if (classNoise > 0) {
        ns<-runif(noInst,0,1)
        Cn <- as.integer(runif(len,1, maxValue+1))
        Cnorm[ns<classNoise]<-Cn # change to random variable
    }
    #output
    out<-as.data.frame(cbind(Cnorm,P[[1]],P[[2]],B[[1]],B[[2]],E[[1]],E[[2]],R[[1]],R[[2]]))
    names(out)<-c("class","Pweak","Pstrong","Bweak","Bstrong","Eweak","Estrong","Iuniform","Inormal")
    out
}

basicStrong<-function(A) {
    X<-vector(mode="integer",length=length(A))
    X[A<=3] <- -4
    X[A==4] <- -2
    X
}
basicWeak<-function(A) {
    X<-vector(mode="integer",length=length(A))
    X[A<=2] <- -2
    X
}
performanceWeak<-function(A){
    X<-vector(mode="integer",length=length(A))
    X[A==1] <- -3
    X[A==2] <- -2
    X[A==4] <- 2
    X[A==5] <- 3
    X   
}
performanceStrong<-function(A){
    X<-vector(mode="integer",length=length(A))
    X[A==1] <- -5
    X[A==2] <- -3
    X[A==4] <- 3
    X[A==5] <- 5
    X   
}
excitementWeak<-function(A){
    X<-vector(mode="integer",length=length(A))
    X[A<=4] <- 0
    X[A==5] <- 1 
    X
}
excitementStrong<-function(A){
    X<-vector(mode="integer",length=length(A))
    X[A<=4] <- 0
    X[A==5] <- 4 
    X
}
distGen<-function(len, probs){
    sumP <- sum(probs)
    normProbs <- probs / sumP
    cumProbs <- c(normProbs[1])
    for (i in 2: length(normProbs))
        cumProbs[i] <- cumProbs[i-1]+normProbs[i]
    rnd <- runif(len,0,1)
    gen <- vector(mode="integer", length=len)
    gen[rnd <= cumProbs[1]] <- 1
    for (i in 2:length(cumProbs)){
        gen[rnd > cumProbs[i-1] & rnd <= cumProbs[i]] <- i
    }
    gen
}

## ties are broken by first come gets lower score
forceDist <- function(A,probs) {
    sumP <- sum(probs)
    normProbs <- probs / sumP
    cumProbs <- c(normProbs[1])
    for (i in 2: length(normProbs))
        cumProbs[i] <- cumProbs[i-1]+normProbs[i]
    cumIdx <- as.integer(cumProbs*length(A))
    cumIdx[length(cumIdx)] <- length(A) ## to avoid errors due to numerical rounding 
    oA<-order(A)
    gen <- vector(mode="integer", length=length(A))
    gen[oA[1:cumIdx[1]]] <- 1
    for (i in 2:length(cumIdx)){
        gen[oA[cumIdx[i-1]:cumIdx[i]]] <- i
    }
    gen    
}
