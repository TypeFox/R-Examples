hmm <- function(y,yval=NULL,par0=NULL,K=NULL,rand.start=NULL,stationary=cis,
                mixture=FALSE,cis=TRUE,tolerance=1e-4,verbose=FALSE,itmax=200,
                crit='PCLL',keep.y=TRUE,data.name=NULL) {
#
# Function hmm.  To fit a Hidden Markov model to a data set where the
# observations come from one of a number of finite discrete
# distributions, depending on the (hidden) state of the Markov chain.
# These distributions are specified by a matrix Rho = [rho_ij] where
# rho_ij = P(X = x_i | S = j), X being the observable random variable
# and S being the hidden state.
# 
# Note that y used to be allowed to be a matrix, each column being interpreted
# as an independent replicate of the observation sequence.  Now y should be
# either a vector or a *list* of vectors.
#

# Check on consistency of ``mixture'' and ``stationary''.
if(mixture & !stationary)
	stop("Makes no sense for mixture to be non-stationary.\n")

# Check on consistancy of ``stationary'' and ``cis''.
if(stationary & !cis)
	stop(paste("If the model is stationary the initial state\n",
                   "probability distribution must be constant\n",
                   "across data sequences.\n"))

# Put together a data name tag for the output.
if(is.null(data.name)) data.name <- deparse(substitute(y))

# Change y into *character* data.  If y is a matrix, change it
# to a list, and put out a snarky message to the user.
y <- charList(y)

# Check that the observation values are compatible
# with yval if it is specified.
uval <- attr(y,"uval")
if(is.null(yval)) yval <- uval
if(!all(uval%in%yval))
        stop("Specified y values do not include all observed y values.\n")
lns  <- sapply(y,length)
nval <- length(yval)

# Check that one of par0 and K is specified and use the specified
# one to determine the other.
if(is.null(par0) & is.null(K))
	stop('One of par0 and K must be specified.')
if(is.null(par0)) {
	if(is.null(rand.start)) rand.start <- list(tpm=FALSE,Rho=FALSE)
	par0  <- init.all(nval,K,rand.start,mixture)
}
else {
	K <- nrow(par0$tpm)
        if(K != ncol(par0$tpm))
            stop("The specified tpm is not square.\n")
	if(nrow(par0$Rho) < nval)
		stop(paste("The row dimension of Rho is less than\n",
                           "the number of distinct y-values.\n"))
}
if(is.null(rownames(par0$Rho))) {
    if(length(yval) != nrow(par0$Rho)) {
        whinge <- paste("No rownames for Rho and nrow(Rho) is not equal\n",
                        "to the number of unique y values.\n",sep="")
        stop(whinge)
    }
        rnms <- rownames(par0$Rho) <- yval
} else {
    rnms <- rownames(par0$Rho)
    if(!all(yval%in%rnms)) {
        whinge <- paste("The row names of the initial value of \"Rho\" do not\n",
                        "include all possible y-values.\n")
        stop(whinge)
    }
}

# If K=1 do the triv thing:
if(K==1) {
	y <- factor(unlist(y),levels=yval)
	Rho <- as.matrix(table(y)/length(y))
        rownames(Rho) <- rnms
	ll  <- sum(log(ffun(y,Rho)))
	return(list(Rho=Rho,tpm=NA,ispd=NA,log.like=ll,
               converged=NA,nstep=NA,data.name=data.name))
}

# Pick out the index of the stopping criterion:
icrit <- match(crit,c('PCLL','L2','Linf'))
if(is.na(icrit)) stop('Stopping criterion not recognized.')

# Perform initial setting-up.
tpm    <- par0$tpm
if(stationary) {
    ispd   <- revise.ispd(tpm)
} else { # Make the chains equally likely to start in any state.
    ispd <- matrix(1/K,K,length(y))
}
Rho    <- par0$Rho
m      <- nrow(Rho)
digits <- 2+ceiling(abs(log10(tolerance)))

old.theta <- c(c(tpm[,-K]),c(Rho[1:(m-1),]))
fy        <- ffun(y,Rho)
rp        <- recurse(fy,tpm,ispd,lns)
old.ll    <- sum(log(rp$llc))

if(verbose) cat('\n      Initial set-up completed ...\n\n')

# Revise:
em.step <- 1
if(verbose) cat('Repeating ...\n\n')
chnge <- numeric(3)
repeat{
	if(verbose) cat(paste('EM step ',em.step,':\n',sep=''))

# Calculate the parameters.
	tpm  <- revise.tpm(rp$xi,mixture)
	ispd <- if(stationary) {
			revise.ispd(tpm)
		} else {
			revise.ispd(gamma=rp$gamma,lns=lns,cis=cis)
		}
	Rho  <- revise.rho(y,rp$gamma,rnms)

# Update the log likelihood on the basis of the
# new parameter estimates.  This entails calculating
# the new recursive probabilities (which will be used
# to update the parameter estimates on the *next* EM
# step, if necessary).
	fy <- ffun(y,Rho)
	rp <- recurse(fy,tpm,ispd,lns)
	ll <-  sum(log(rp$llc))

# Test for convergence:
	new.theta <- c(c(tpm[,-K]),c(Rho[1:(m-1),]))
	chnge[1] <- 100*(ll - old.ll)/abs(old.ll)
	chnge[2] <- sqrt(sum((old.theta-new.theta)**2))
	chnge[3] <- max(abs(old.theta-new.theta))
	if(verbose){
		cat('     Log-likelihood: ',
            	format(round(ll,digits)),'\n',sep='')
		cat('     Percent decrease in log-likelihood: ',
		    format(round(chnge[1],digits)),'\n',sep='')
		cat('     Root-SS of change in coef.: ',
		    format(round(chnge[2],digits)),'\n',sep='')
		cat('     Max. abs. change in coef.: ',
		    format(round(chnge[3],digits)),'\n',sep='')
	}

	if(chnge[icrit] < tolerance) {
			converged <- TRUE
			nstep <- em.step
			break
		}

	if(em.step >= itmax) {
		cat('Failed to converge in ',itmax,' EM steps.\n',sep='')
		converged <- FALSE
		nstep <- em.step
		break
	}

# Replace the ``old'' parameter and log likelihood values
# by the new ones.
	old.theta <- new.theta
	old.ll    <- ll
# Increment the step number.
	em.step   <- em.step + 1
}

# Tidy up a bit:
if(length(y)==1) {
   if(keep.y) y <- y[[1]]
   ispd <- as.vector(ispd)
}
stnms <- rownames(par0$tpm)
if(!is.null(stnms)) {
    rownames(tpm) <- stnms
    colnames(tpm) <- stnms
    names(ispd)   <- stnms
    colnames(Rho) <- stnms
}
rslt <- list(Rho=Rho,tpm=tpm,ispd=ispd,log.like=ll,converged=converged,
             nstep=nstep,y=if(keep.y) y else NULL, data.name=data.name,
             stationary=stationary)
class(rslt) <- "hmm.discnp"
rslt
}
