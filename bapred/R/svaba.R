svaba <-
function(x, y, batch, nbf=NULL, algorithm="fast") {

  if(any(is.na(x)))
	stop("Data contains missing values.")
  if(!is.matrix(x))
    stop("'x' has to be of class 'matrix'.") 
  if(!(is.factor(y) & all(levels(y)==(1:2))))
    stop("'y' has to be of class 'factor' with levels '1' and '2'.")

  ##require("sva")

  if(!is.null(nbf)) {
    numsv <- nbf
  }
  else {
  # There has to be a case differentiation, because the estimation of the
  # number of factors using 'num.sv' does only work, when the number of
  # variables is at least the number of observations:
  if(nrow(x) <= ncol(x))
      numsv <- sva::num.sv(dat=t(x), mod=cbind(1, as.numeric(y)-1))
    else
      numsv <- sva::num.sv(dat=t(x[1:ncol(x),]), mod=cbind(1, as.numeric(y[1:ncol(x)])-1))
  }
  
  if (numsv!=0) {
	
    svobj <- sva::sva(dat = t(x), mod = cbind(1, as.numeric(y) - 
        1), n.sv = numsv)
			
    mod <- cbind(1, as.numeric(y) - 1)
    nmod <- dim(mod)[2]
    mod <- cbind(mod, svobj$sv)
    gammahat <- (t(x) %*% mod %*% solve(t(mod) %*% mod))[, (nmod + 
        1):(nmod + numsv)]
			
    db = t(x) - gammahat %*% t(svobj$sv)
    xadj <- t(db)
		
    params <- list(xadj = xadj, xtrain = x, ytrain = y, svobj = svobj, 
        algorithm = algorithm)
		
  }
  else {

    warning("Estimated number of factors was zero.")  
    xadj <- x
    params <- list(xadj = xadj, xtrain = x, ytrain = y, svobj = NULL, 
        algorithm = algorithm)
		
  }    
	
  params$nbatches <- length(unique(batch))
  params$batch <- batch
	
  class(params) <- "svatrain"

  return(params)

}
