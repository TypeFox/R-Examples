class.Lee <-
function( cutscore, ip, ability = NULL, rdm = NULL, quadrature=NULL, D = 1.7){

		if(is.null(quadrature)){
		if(is.null(ability)){

			theta <- MLE(rdm, ip, D)} else{
			theta <- ability}

		results <- Lee.P(cutscore, ip, theta, D)
		results

		} else {

			if(length(quadrature)!=2) stop("quadrature points and weights must be a list of length 2")
			if(length(quadrature[[1]])!=length(quadrature[[2]])) stop("number of quadrature points and weights do not match")

		results <- Lee.D(cutscore, ip, quadrature, D)

	results}
  }
