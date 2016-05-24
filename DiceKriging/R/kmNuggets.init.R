`kmNuggets.init` <-
function(model) {

	n <- nrow(model@X)
	parinit <- model@parinit
	lower <- model@lower
	upper <- model@upper
	ninit <- model@control$pop.size
	param.n <- model@covariance@param.n
  
	if (model@covariance@nugget.flag & !model@covariance@nugget.estim) {
    nugget.aux <- rep(model@covariance@nugget, n)
	}
  
  if (length(parinit)>0) {
	  matrixinit <- matrix(parinit, nrow = param.n, ncol = ninit) 
	} else {
	  if (existsMethod("paramSample", signature = class(model@covariance))) {
      matrixinit <- paramSample(model@covariance, n=ninit, lower=lower, upper=upper, y=model@y)
    } else {
      # sample ninit design points, generated from uniform [lower, upper]
      matrixinit <- matrix(runif(ninit*param.n), nrow = param.n, ncol = ninit)
      matrixinit <- lower + matrixinit*(upper - lower)
	  }
	}
  
  if (model@noise.flag) nugget.aux <- model@noise.var
  
	
		# variance standard estimate (biased negatively)
	trend.estimate <- lm(model@y~model@F-1)
	random.part.estimate <- trend.estimate$residuals
	varinit.total <- var(random.part.estimate)
	varinit.standard <- varinit.total - mean(nugget.aux)
	if (varinit.standard <=1e-20) varinit.standard <- 1/2*varinit.total
	
		# variance estimate using the variogram (biased negatively, but less)
	x.dist <- dist(model@X)
	y.dist <- dist(random.part.estimate)
	I <- (x.dist > quantile(x.dist, 0.5))
	matrix.nugget.aux <- matrix(nugget.aux, n, n)
	matrix.sym.nugget.aux <- (matrix.nugget.aux + t(matrix.nugget.aux))/2
	nugget.aux.I <- matrix.sym.nugget.aux[I]
	varinit.vario.total <- 1/2*mean(y.dist[I]^2)
	varinit.vario <- varinit.vario.total - mean(nugget.aux.I)
	if (varinit.vario<=1e-20) varinit.vario <- 1/2*varinit.vario.total
	
		# final choice
	varinit <- (varinit.standard + varinit.vario)/2
	
	
		# boundaries    
	varinit.upper <- varinit.total - min(nugget.aux) 
	varinit.lower <- varinit.total - max(nugget.aux)
	if (varinit.upper<=1e-20) varinit.upper <- varinit.total
	if (varinit.lower<=1e-20) varinit.lower <- 1e-20
			
	varinit.vario.upper <- varinit.vario.total - min(nugget.aux.I) 
	varinit.vario.lower <- varinit.vario.total - max(nugget.aux.I)
	if (varinit.vario.upper<=1e-20) varinit.vario.upper <- varinit.vario.total
	if (varinit.vario.lower<=1e-20) varinit.vario.lower <- 1e-20
	
		# final choice
	varinit.lower <- 1/10*min(varinit.lower, varinit.vario.lower)
	varinit.upper <- 10*max(varinit.upper, varinit.vario.upper)
			
	varinit.sim <- runif(n=ninit, min=1/2*varinit, max=3/2*varinit)       
					
    # take the best point(s)				 
	matrixinit <- rbind(matrixinit, varinit.sim)
	fninit <- apply(matrixinit, 2, logLikFun, model)
	selection <- sort(fninit, decreasing = TRUE, index.return = TRUE)$ix
	selection <- selection[1:model@control$multistart]
	parinit <- matrixinit[, selection, drop = FALSE]
  # for one point : parinit <- matrixinit[, which.max(logLikinit)] 
	
	lp <- nrow(parinit)
	covinit <- list()
	for (i in 1:model@control$multistart){
    pari <- as.numeric(parinit[, i])
	  covinit[[i]] <- vect2covparam(model@covariance, pari[1:(lp-1)])
	  covinit[[i]]@sd2 <- pari[lp]
	}
	
	return(list(par = parinit, 
	            value = fninit[selection], 
	            cov = covinit,
	            lower = c(lower, varinit.lower), 
	            upper = c(upper, varinit.upper)))
  
}

