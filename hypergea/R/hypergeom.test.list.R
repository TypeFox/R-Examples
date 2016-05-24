hypergeom.test.list <-
function( x, alternative='two.sided', pval.method="fisheR", nthreads=2, ... ){
	getProbObs <- function(margins, CT, lfN, len_dm){
 		lm0 <- .Call("lfactorial", (c(CT)), PACKAGE="hypergea" )
 		lfm <- .Call("lfactorial", ((margins)), PACKAGE="hypergea" )
		
		dens <- exp( lfm - ( (len_dm-1)*lfN + lm0 ) )
	return(dens)
	}

	x <- lapply(x, as.array)
	dm <- as.integer(dim(x[[1]]))
	len_dm <- length(dm)
	CT <- x
	len_CT <- length(CT)
	#CT <- lapply(CT, function(x){apply(x, 1:len_dm, as.integer)})
	
	if( dm[1] == 2 && dm[1]<dm[2] ){ CT <- lapply(CT, t) }
	if( dm[1]<dm[2] ){ CT <- lapply(CT, t) }
	
	N <- as.integer(sum(CT[[1]])); lfN <- sum(log(1:N))
	dm <- as.integer(dim(CT[[1]]))
	margins <- lapply( 1:len_CT, function(i){ (unlist(getMargins(CT[[i]]))) } )
	prob.obs <- sapply( 1:len_CT, function(i){ getProbObs( margins[[i]], CT[[i]], lfN, len_dm ) } )
	O000 <- vapply( 1:len_CT, function(i){  CT[[i]][1] }, FUN.VALUE=0L )
 	nthreads <- as.integer(nthreads)

	res <- .Call("hypergeom_IxJ_list", O000, N, margins, prob.obs, dm, nthreads, PACKAGE="hypergea")
	res <- lapply( 1:len_CT, function(i){ 
		aux <- c("sum", "less", "greater", "two.sided_ml", "two.sided_fisheR", "two.sided_d")
		res[[i]] <- setNames(as.list(res[[i]]), c("n0", "Prob", "Freq"))
		attr(res[[i]][['Prob']], "names") <- aux
		attr(res[[i]][['Freq']], "names") <- aux
		return(res[[i]])
	} )
	or <- rep(NA, times=len_CT); if(dm[1]==2 && dm[2]==2){or <- sapply(CT, getOddsRatio)}

 	cse <- sapply(1:len_CT, function(i){ 1.96*sqrt( sum(1/(CT[[i]]+ifelse(any(CT[[i]]==0), 0.5, 0))) ) })
 	conf.int <- cbind(exp(log(or)-cse), exp(log(or)+cse))
 	conf.int <- lapply( 1:len_CT,  function(i){ x <- conf.int[i,]; attr(x, "conf.level") <- 0.95; return(x) })

	p.value <- rep(1, times=len_CT);
	if(alternative %in% c("less", "greater")){
		p.value <- sapply(res, function(x){x$Prob[alternative]})
		pval.method <- "exact"
	}else{
		pv <- sapply(1:len_CT, function(i){
			x <- res[[i]]$Prob
			p.value0 <- switch(pval.method, fisheR=x['two.sided_fisheR'], minimum.likelihood = x['two.sided_ml'], 
					double = x['two.sided_d'])
			if( p.value0 == 0.0 && pval.method != 'two.sided_fisheR' ){ # p.value equal to zero can occur when, e.g. x[1,0]==0 and x[0,0] small (<<10)
				p.value0 <- x['two.sided_ml']
				if( p.value0 == 0.0 ){ 
					p.value0 <- x['two.sided_d'] 
				} 
			}
			return(p.value0)
		} )
		p.value <- pv
	}
	
	METHOD <- "Exact hypergeometric test for a "
	if( length(dm) == 3 ){ METHOD <- paste(METHOD, "2x2x2 table", sep="") }
	if( length(dm) == 2 ){ I <- dm[1]; J <- dm[2]; METHOD <- paste(METHOD, I, "x", J, " table", sep="") }

	OBJ <- lapply( 1:len_CT, function(i){
		obj <- list(p.value=as.numeric(p.value[i]), pval.method=as.character(pval.method)
		, alternative=alternative, statistic=setNames(O000[i], "a"), estimate=setNames(or[i], "odds ratio"), prob.obs=prob.obs[i], method=METHOD
		, data.name=paste(c(deparse(substitute(x)), res$call), collapse=", "), conf.int=conf.int[[i]])
		class(obj) <- "htest"
		return(obj)
	})
	names(OBJ) <- names(CT)

return(OBJ)
}
