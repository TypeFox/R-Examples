

NMIextract <- function(results, expr, fun){
  pf<-parent.frame()
  results0 <- results
  if (!is.null(match.call()$expr)){
    expr<-substitute(expr)
	lapply( results0 , FUN = function(results){
        lapply(results, function(result) eval(expr, result,pf))
							} )
  } else {
    lapply( results0 , FUN = function( results){
				lapply(results, fun)
							} )
  }
  
}