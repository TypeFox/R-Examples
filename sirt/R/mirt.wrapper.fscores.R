
##################################################################
# calculation of factor scores (EAP, MAP and MLE) in mirt
# Note that MAP and MLE are evaluated at the discrete grid
mirt.wrapper.fscores <- function( mirt.obj , weights=NULL ){
	# posterior distribution of mirt object
	mirt.post <- mirt.wrapper.posterior(mirt.obj , weights=weights)	
	data <- mirt.post$data
	D <- ncol(mirt.post$theta.k)
	theta.k <- mirt.post$theta.k
	if ( is.null(weights) ){ 
		N <- nrow(data) 
		weights <- rep(1,N)
				}
	p.xi.aj <- mirt.post$f.yi.qk
	p.aj.xi <- mirt.post$f.qk.yi
	res <- .smirt.person.parameters( data=data , D=D , theta.k=theta.k ,
		   p.xi.aj , p.aj.xi , weights=weights )
	return(res)
	}
##################################################################	