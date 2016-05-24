## ============================
## Matrix of substitution costs
## ============================

seqsubm <- function(seqdata, method, cval=NULL, with.missing=FALSE, 
	miss.cost=NULL, time.varying=FALSE, weighted=TRUE, transition="both", lag=1, missing.trate=FALSE) {

	if (!inherits(seqdata,"stslist"))
		stop(" [!] data is NOT a sequence object, see seqdef function to create one")

	metlist <- c("CONSTANT","TRATE")
	if (missing(method) || !method %in% metlist)
		stop(" [!] method must be one of: ", paste(metlist,collapse=" "))
	
	transitionlist <- c("previous", "next", "both")
	if (!transition %in% transitionlist)
		stop(" [!] transition must be one of: ", paste(transitionlist,collapse=" "))

	alphabet <- attr(seqdata,"alphabet")
	
	cval4cond <- time.varying && method=="TRATE" && transition=="both"
	if (is.null(cval)) {
		cval <- ifelse(cval4cond, 4, 2)
	}
	if (is.null(miss.cost)) {
		miss.cost <- cval
	}
	## Adding an entry for for missing state
	if (with.missing) {
		if(missing.trate && method == "TRATE") {
			message(" [>] using transition rate to compute substitution cost for missing values")
		} else {
			message(" [>] setting ",miss.cost," as substitution cost for missing values")
		}
		
		alphabet <- c(alphabet,attr(seqdata,"nr"))
	}

	alphsize <- length(alphabet)
	
	if (method=="CONSTANT") {
		if (is.na(cval)){
			stop("no value for the constant substitution-cost")
		}
		if (time.varying) {
			time <- ncol(seqdata)
			
			message(" [>] creating ",alphsize,"x",alphsize,"x",time,
				" time varying substitution-cost matrix using ",
				cval," as constant value")
			costs <- array(cval, dim=c(alphsize, alphsize, time))
			for (i in 1:time) {
				diag(costs[,,i]) <- 0
			}
		}
		else {
			message(" [>] creating ",alphsize,"x",alphsize,
				" substitution-cost matrix using ",cval," as constant value")

			costs <- matrix(cval,nrow=alphsize,ncol=alphsize)
			diag(costs) <- 0
		}
	}

	if (method=="TRATE") {
		if (time.varying) {
			message(" [>] creating time varying substitution-cost matrix using transition rates ...")
			tr <- seqtrate(seqdata, time.varying=TRUE, weighted=weighted, lag=lag, with.missing=missing.trate)
			tmat <- nrow(tr)
			time <- ncol(seqdata)
			costs <- array(0, dim=c(alphsize, alphsize, time))
			## Function to compute the cost according to transition rates
			tratecostBoth <- function(trate, time, state1, state2, debut, fin){
				cost <- 0
				if (!debut) { ## Premier état
					cost <- cost - trate[state1,state2, time-1] - trate[state2, state1, time-1]
				}
				if (!fin) { ##Dernier Etat
					cost <- cost - trate[state1,state2, time] - trate[state2, state1, time]
				}
				if (!debut && !fin) {
					return(cost + cval)
				}
				else{
					return(cval + 2*cost)
				}
			}
			tratecostPrevious <- function(trate, time, state1, state2, debut, fin){
				cost <- 0
				if (!debut) { ## Premier état
					cost <- cost - trate[state1,state2, time-1] - trate[state2, state1, time-1]
				}
				return(cval + cost)
			}
			tratecostNext <- function(trate, time, state1, state2, debut, fin){
				cost <- 0
				if (!fin) { ##Dernier Etat
					cost <- cost - trate[state1,state2, time] - trate[state2, state1, time]
				}
				return(cval + cost)
			}
			if(transition=="previous")	{
				tratecost <- tratecostPrevious
			}
			else if (transition=="next") {
				tratecost <- tratecostNext
			}
			else {
				tratecost <- tratecostBoth
			}
			
			for (t in 1:time) {
				for (i in 1:(tmat-1)) {
					for (j in (i+1):tmat) {
						cost <- max(0,tratecost(tr, t, i, j, t==1, t==time))
						costs[i, j, t] <- cost
						costs[j, i, t] <- cost
					}
				}
			}
		}
		else {
			message(" [>] creating substitution-cost matrix using transition rates ...")
			tr <- seqtrate(seqdata, time.varying=FALSE, weighted=weighted, lag=lag, with.missing = missing.trate)
			tmat <- nrow(tr)
			costs <- matrix(nrow=alphsize,ncol=alphsize)
			diag(costs) <- 0
			for (i in 1:(tmat-1)) {
				for (j in (i+1):tmat) {
					cost <- cval - tr[i,j] - tr[j,i]
					costs[i,j] <- cost
					costs[j,i] <- cost
				}
			}
		}
	}

	## 
	if (with.missing &&(method=="CONSTANT"||(method=="TRATE" && missing.trate==FALSE))) {
		if (time.varying) {
			costs[alphsize,1:(alphsize-1),] <- miss.cost
			costs[1:(alphsize-1),alphsize,] <- miss.cost
		}
		else {
			costs[alphsize,1:(alphsize-1)] <- miss.cost
			costs[1:(alphsize-1),alphsize] <- miss.cost
		}
	}

	## Setting rows and columns labels
	rclab <- paste(alphabet,"->",sep="")
	if (time.varying) {
		dimnames(costs) <- list(rclab, rclab, colnames(seqdata))
	}
	else {
		dimnames(costs) <- list(rclab,rclab)
	}

	return(costs)
}

