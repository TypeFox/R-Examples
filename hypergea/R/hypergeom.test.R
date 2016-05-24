hypergeom.test <-
function( x, alternative='two.sided', pval.method="fisheR", nthreads=2, ... ){
	res <- .hypergeax.test( x=x, nthreads=nthreads, ... )
	dm <- dim(x)

	p.value <- 1;
	if(alternative %in% c("less", "greater")){
		p.value <- res$prob[alternative]
		pval.method <- "exact"
	}else{
		p.value <- switch(pval.method, fisheR=res$prob['two.sided_fisheR'], minimum.likelihood = res$prob['two.sided_ml'], double = res$prob['two.sided_d'])
		if( p.value == 0.0 && pval.method != 'two.sided_fisheR' ){ # p.value equal to zero can occur when, e.g. x[1,0]==0 and x[0,0] small (<<10)
			p.value <- res$prob['two.sided_ml']
			if( p.value == 0.0 ){ 
				p.value <- res$prob['two.sided_d'] 
			} 
		} 
	}

	
	METHOD <- "Exact hypergeometric test for a "
	if( length(dm) == 3 ){ METHOD <- paste(METHOD, "2x2x2 table", sep="") }
	if( length(dm) == 2 ){ I <- dm[1]; J <- dm[2]; METHOD <- paste(METHOD, I, "x", J, " table", sep="") }

	obj <- list(p.value=as.numeric(p.value), pval.method=as.character(pval.method), alternative=alternative, statistic=setNames(res$count.obs, "a"), estimate=setNames(res$odds.ratio, "odds ratio"), prob.obs=res$prob.obs, method=METHOD
	, data.name=paste(c(deparse(substitute(x)), res$call), collapse=", "), conf.int=res$conf.int)
	class(obj) <- "htest"

	if( !is.null(list(...)$r) ){ obj <- list(obj, res) }
return(obj)
}
