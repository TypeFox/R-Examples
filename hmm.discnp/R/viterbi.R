viterbi <- function(y,object=NULL,tpm,Rho,ispd=NULL,log=FALSE) {
#
# Function viterbi to apply the Viterbi algorithm to a collection
# of data sequences, given the parameters of the model.
#

# If ``object'' is present, get the parameters from that, and
# ignore those (if any) specified as separate arguments.
if(!is.null(object)) {
	tpm  <- object$tpm
	Rho  <- object$Rho
	ispd <- object$ispd
}
K <- nrow(tpm)
if(missing(y)) {
	y <- if(!is.null(object)) object$y else NULL
	if(is.null(y)) stop("No observation sequence supplied.\n")
}
y <- charList(y)

# Build ispd if it was given as NULL
if(is.null(ispd)) ispd <- revise.ispd(tpm)

# Make sure that the y-values are compatible with Rho.
Rho <- check.yval(y,Rho)

# Make sure y is a list, and get the number of sequences and
# lengths of these sequences.
if(is.atomic(y)) y <- list(y)
nseq <- length(y)
lns  <- sapply(y,length)

rslt <- list()
for(j in 1:nseq) {
	psi <- list()
        if(log) {
            delta <- log(ispd) + log(Rho[y[[j]][1],])
        } else {
	    delta <- ispd*Rho[y[[j]][1],]
	    delta <- delta/sum(delta)
        }
	nj <- lns[j]
	for(tt in 2:nj) {
		if(log) {
                    tmp <- apply(delta + log(tpm),2,
                             function(x){((1:length(x))[x==max(x)])}
                             )
                } else {
		    tmp <- apply(delta*tpm,2,
                             function(x){((1:length(x))[x==max(x)])}
                             )
                }
	        psi[[tt]] <- tmp # Note that tmp will be a list of
		                 # vectors, each of length between
                                 # 1 and K = the number of states.
		if(log) {
                    delta <- log(Rho[y[[j]][tt],]) +
                                 apply(delta + log(tpm),2,max)
                } else {
		    delta <- Rho[y[[j]][tt],]*apply(delta*tpm,2,max)
                    delta <- delta/sum(delta)
                }
	}
	temp <- list()
	temp[[nj]] <- (1:K)[delta==max(delta)]
	for(tt in (nj-1):1) {
		i <- 0
		temp[[tt]] <- list()
		for(x in temp[[tt+1]]) {
			k <- x[1]
			for(w in psi[[tt+1]][[k]]) {
				i <- i+1
				temp[[tt]][[i]] <- c(w,x)
			}
		}
	}
        rrr <- matrix(unlist(temp[[1]]), nrow = nj)
        rslt[[j]] <- if(ncol(rrr)==1) as.vector(rrr) else rrr
}
if(nseq==1) rslt[[1]] else rslt
}
