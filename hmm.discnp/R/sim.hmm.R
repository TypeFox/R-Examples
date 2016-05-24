sim.hmm <- function(nsim,tpm,Rho,ispd=NULL,yval=NULL,verb=FALSE) {
#
# Function sim.hmm to simulate data from a hidden Markov
# model with transition probability matrix tpm, and discrete
# (non-parametric) distributions specified by the matrix Rho.
#

# Check for validity of the Rho argument:
if(any(Rho<0)) stop("Negative entries in Rho.\n")
Rho <- as.matrix(Rho)
xxx <- unname(apply(Rho,2,sum))
if(!isTRUE(all.equal(xxx,rep(1,ncol(Rho)))))
	stop("Columns of Rho do not all sum to 1.\n")

# Check for validity of yval argument:
if(is.null(yval)) {
	yval <- 1:nrow(Rho)
} else {
	if(length(yval)!=nrow(Rho))
		stop(paste("Mismatch between length of \"yval\"",
                           "and number of rows of \"Rho\".\n"))
}
row.names(Rho) <- yval
nseq <- length(nsim)

# If Rho has a single column, generate i.i.d. data.
if(ncol(Rho)==1) {
	temp <- lapply(nsim,function(n,yval,Rho){
				sample(yval,size=n,prob=Rho[,1],replace=TRUE)
                            },yval=yval,Rho=Rho)
	return(if(nseq==1) temp[[1]] else temp)
}

if(ncol(tpm) != nrow(tpm))
	stop("The matrix tpm must be square.\n")
if(ncol(tpm) != ncol(Rho))
	stop("Mismatch between dimensions of tpm and Rho.\n")
if(any(tpm<0)) stop("Negative entries in tpm.\n")
xxx <- unname(apply(tpm,1,sum))
if(!isTRUE(all.equal(xxx,rep(1,nrow(tpm)))))
	stop("Rows of tpm do not all sum to 1.\n")

if(is.null(ispd)) ispd <- revise.ispd(tpm)
K    <- ncol(Rho)
M    <- nrow(Rho)
ntot <- sum(nsim)
rslt <- if(is.numeric(yval)) numeric(ntot) else character(ntot)
jr   <- 0
for(j in 1:nseq) {
	jr <- jr+1
	s1     <- sample(1:K,1,prob=ispd)
        rslt[jr] <- sample(yval,1,prob=Rho[,s1])
	for(i in 2:nsim[j]) {
		jr <- jr+1
		s1 <- sample(1:K,1,prob=tpm[s1,])
		rslt[jr] <- sample(yval,1,prob=Rho[,s1])
		if(verb) {
			if(jr%%1000 == 0) cat(jr,"")
			if(jr%%10000 == 0) cat("\n")
		}
	}
}
if(verb) cat("\n")
if(nseq==1) rslt else {
	rslt <- unname(tapply(rslt,rep(1:nseq,nsim),function(x){x}))
	attr(rslt,"dim") <- NULL
	rslt
}
}
