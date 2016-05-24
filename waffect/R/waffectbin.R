waffectbin = function(prob, count, label, method, burnin){
	r = as.integer(count[1]) 
	ninds = length(prob)
	#Call C++ function waffectbin
  	if (method=="mcmc") {
          res <- .Call( "waffectbin_mcmc", prob , r , as.integer(burnin) , PACKAGE = "waffect" )
        } else if (method=="reject") {
          res <- .Call( "waffectbin_reject", prob , r , PACKAGE = "waffect" )
        } else {
          res <- .Call( "waffectbin", prob , r , as.integer(ninds), PACKAGE = "waffect" )
        }

	# Affect the labels
	return(label[(!res)+1]);
}




