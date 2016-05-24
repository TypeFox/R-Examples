io <-
function(...,data=NULL) 
{
  # This function can take a bunch of numeric sample vectors 
  # in ..., or a list of such, or a formula that specifies a
  # a response (e.g., y), grouped by a factor, e.g., g, of same
  # length as y, via y ~ g. 
  # This breaks down y into the desired samples, one sample 
  # for each factor level. 
    xlist <- list(...)
    if(is(xlist[[1L]], "formula")) {
	cl <- match.call()  # gets a copy of the current call
    	mf <- cl 
    	mf[[1L]] <- as.name("model.frame") 
    	mf <- eval(mf,parent.frame())
    # mf is a data frame consisting of response 
    # and the other model variable
    	mt <- attr(mf, "terms")
    	y <- model.response(mf, "numeric") # response vector
    	fname <- attr(mt,"term.labels") 
    # fname contains the names of explanatory variables
    	if(length(fname) != 1) {
      	stop("Please specify only one term in the formula")
    	}
    	fvec <- as.factor(mf[, fname])
    # fvec contains the values of the single explanatory variable
    # as a factor
    samples <- lapply(levels(fvec), 
                      function(flvl) {
                        return(y[which(fvec == flvl)])} )
  } else {
    # tests whether ... is a list, when not a formula
    	if (is.list(xlist[[1]])) {
      	samples <- xlist[[1]]
    	} else {
	if( all(unlist(lapply(xlist,FUN=is.numeric))) == FALSE) stop("improper input for ...\n")
        samples <- lapply(xlist,as.numeric)
     }
    
  }
  if(length(samples) < 2) stop("fewer than 2 samples\n")
  if(all(unlist(lapply(samples,FUN=is.numeric))) != TRUE) stop("improper input for ...\n")
  lapply(samples,as.numeric)
  
}
