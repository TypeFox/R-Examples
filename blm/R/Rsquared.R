Rsquared <- function(object){

if(class(object)[1]!="lexpit"&class(object)[1]!="blm")
		stop("Object must be an instance of a blm or lexpit model.")
		
mcfadden <- function(loglik,loglik.null){
  1-loglik/loglik.null
}

mcfadden.adj <- function(loglik,loglik.null,num.params){
  1-(loglik-num.params)/loglik.null
}
	
	if(class(object)=="lexpit")
		num <- object@p+object@q
	else
		num <- length(object@coef)
		
	list(
		R2 = mcfadden(object@loglik,object@loglik.null),
		R2adj = mcfadden.adj(object@loglik,object@loglik.null,num)
	 )
	
}