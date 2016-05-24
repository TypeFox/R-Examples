.hypergeax.test <-
function( x, nthreads=2, ... ){
	getProbObs <- function(margins, CT, N, dm){
 		lm0 <- .Call("lfactorial", c(CT), PACKAGE="hypergea" )
 		lfm <- .Call("lfactorial", unlist(margins), PACKAGE="hypergea" )

 		lfN <- sum(log(1:N))
		dens <- exp( lfm - ( (length(dm)-1)*lfN + lm0 ) )
	return(dens)
	}

	x <- as.array(x)
	dm <- as.integer(dim(x))
	len_dm <- length(dm)

	if( any(is.na(x)) ){ stop(paste0("Unable to deal with NA in ", sQuote("x"))) }
	CT <- apply(x, 1:len_dm, as.integer)
	if( !all( CT == x ) ){ warning(paste0(sQuote("x"), " contains non-integer values, are truncated")) }
	if( !(len_dm < 4) ){ warning(paste0("dimension of ", sQuote("x"), " currently not supported")) }
	
	if( dm[1] == 2 && dm[1]<dm[2] ){ CT <- t(CT) }
	if( dm[1]<dm[2] ){ CT <- t(CT) }
	
	N <- as.integer(sum(CT))
	dm <- as.integer(dim(CT))
	margins <- getMargins( CT )
	prob.obs <- getProbObs( margins, CT, N, dm )
	marg.obs <- as.integer(sapply(margins, function(x){x[1]}))
	Prob <- double(6)
	Freq <- double(6)
	n0 <- 0
	or <- 0
	Hist <- double(0)

 	nthreads <- as.integer(nthreads)

	res <- c()
	if( (len_dm == 2) ){
		if( !is.null(list(...)$aa) ){ 		
		res <- .Call("hypergeom_IxJ_a", CT[1], N, (unlist(margins)), (prob.obs), dm, nthreads, PACKAGE="hypergea")
		names(res) <- c("n0", "Prob", "Freq")

		}else{
			res <- .C("hypergeom_IxJ", CT[1], N, as.integer(unlist(margins)), as.numeric(prob.obs), n0=as.numeric(n0), Prob=Prob, Freq=Freq, dm, nthreads, PACKAGE="hypergea" )

		}
		or <- NA; if(dm[1]==2 && dm[2]==2){or <- getOddsRatio(CT)}
	}
	if( (len_dm == 3) && all(dm == 2) ){
		marg.obs <- sort(marg.obs, decreasing=FALSE)
		res <- .C("hypergeom_2x2x2", CT[1], N, marg.obs, as.integer(unlist(margins)), as.numeric(prob.obs), n0=as.numeric(n0), Prob=Prob, Freq=Freq, nthreads, PACKAGE="hypergea" )
		or <- getOddsRatio(CT)
	}
	names(res[[ 'Prob' ]]) <- names(res[[ 'Freq' ]]) <- c("sum", "less", "greater", "two.sided_ml", "two.sided_fisheR", "two.sided_d")

	cse <- 1.96*sqrt( sum(1/(CT+ifelse(any(CT==0), 0.5, 0))) )
	conf.int <- exp(c(log(or)-cse, log(or)+cse))
	attr(conf.int, "conf.level") <- 0.95

	obj <- list(
		table=CT, margins=margins, prob=res[[ 'Prob' ]], count=res[[ 'Freq' ]],
		prob.obs=as.numeric(prob.obs), count.obs=as.integer(CT[1]), freq.count.obs=as.numeric(res[[ 'n0' ]]), odds.ratio=or
		, nthreads=nthreads, conf.int=conf.int, hist=Hist
	)
	#if( abs(res[[ 'Prob' ]][1] - 1) > 1e-6 ){ warning( paste("Invalid total probability: ", res[[ 'Prob' ]][1], sep="") ) }

return(obj)
}
