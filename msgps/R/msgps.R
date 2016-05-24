msgps <-
function(X,y,penalty="enet", alpha=0, gamma=1, lambda=0.001,  tau2, STEP=20000, STEP.max=200000, DFtype="MODIFIED", p.max=300, intercept=TRUE, stand.coef=FALSE){

	  #check X,y
	   if(!is.matrix(X)) stop(" X must be a matrix.")
	   if(mode(X)!="numeric") stop(" X must be numeric.")
		if (!is.vector(y)) stop("y must be a vector.")
	   if(mode(y)!="numeric") stop(" y must be numeric.")
		if (nrow(X)!=length(y)) stop("The number of sample must not be differenet between X and y.")
		if (sum(complete.cases(X)==FALSE)>0)  stop("X must be complete data.")
		if (sum(complete.cases(y)==FALSE)>0)  stop("y must be complete data.")

	  #check penalty
		candidate <- c("enet","genet","alasso")
	 	if(sum(candidate == penalty) != 1)	stop('penalty must equal "enet", "genet" or "alasso".')

	  #check alpha, gamma, lambda
	    penalty_int <- which(candidate == penalty)
		penalty_int <- as.integer(penalty_int - 1)
		if(!is.numeric(alpha))	stop('"alpha" must be a numeric.')
		if(!is.numeric(gamma))	stop('"gamma" must be a numeric.')
		if(!is.numeric(lambda))	stop('"lambda" must be a numeric.')
		if(!is.vector(alpha))	stop('"alpha" must be a scalar (1-dimensional vector).')
		if(!is.vector(gamma))	stop('"gamma" must be a scalar (1-dimensional vector).')
		if(!is.vector(lambda))	stop('"lambda" must be a scalar (1-dimensional vector).')
		if(length(alpha) > 1)	stop('"alpha" must be a scalar (1-dimensional vector).')
		if(length(gamma) > 1)	stop('"gamma" must be a scalar (1-dimensional vector).')
		if(length(lambda) > 1)	stop('"lambda" must be a scalar (1-dimensional vector).')
		if(penalty_int==0  &&   (alpha < 0 || alpha >=1))	stop('"alpha" must be in [0,1).')
		if(penalty_int==1  &&   (alpha <= 0 || alpha >=1))	stop('"alpha" must be in (0,1).')
		if(penalty_int==2  &&   gamma<=0 ) stop("gamma must be positive.")
		if(penalty_int==2  &&   lambda<0 ) stop("lambda must be non-negative.")
		if(penalty_int==2  &&   nrow(X) <= ncol(X) && lambda<=0 ) stop("lambda must be non-negative.")

	  #check tau2
		if(missing(tau2)==FALSE){
			if(!is.numeric(tau2))	stop('"tau2" must be a numeric.')
		    if(length(tau2) > 1)	stop('"tau2" must be a scalar (1-dimensional vector).')
			if(tau2<0)	stop('"tau2" must be positive.')
			if(tau2>var(y)*(length(y)-1))	stop('"tau2" must not exceed the variance of y.')
		}

	  #check STEP
		if(mode(STEP)!="numeric")	stop('"STEP" must be numeric.')
		if(length(STEP) > 1)	stop('"STEP" must be a scalar (1-dimensional vector).')
		if(as.integer(STEP)!=STEP)	stop('"STEP" must be integer.')
		if(STEP < 500)	stop('"STEP" must be greater than or equal to 500.')
		if(STEP >= 1e+8)	stop('"STEP" must be less than 1e+8.')

	  #check STEP.max
		if(mode(STEP.max)!="numeric")	stop('"STEP.max" must be numeric.')
		if(length(STEP.max) > 1)	stop('"STEP.max" must be a scalar (1-dimensional vector).')
		if(as.integer(STEP.max)!=STEP.max)	stop('"STEP.max" must be integer.')
		if(mode(STEP.max)!="numeric")	stop('"STEP.max" must be numeric.')
		if(STEP.max< 500)	stop('"STEP.max" must be greater than or equal to 500.')
		if(STEP.max >= 1e+8)	stop('"STEP.max" must be less than 1e+8.')
		if(STEP >= STEP.max)	stop('"STEP.max" must be greater than "STEP".')


	  #check DFtype
		candidate_DFtype <- c("NAIVE","MODIFIED")
	 	if(sum(candidate_DFtype == DFtype) != 1)	stop('DFtype must be "MODIFIED"  or "NAIVE".')

	  #check p.max
		if(mode(p.max)!="numeric")	stop('"p.max" must be numeric.')
		if(length(p.max) > 1)	stop('"p.max" must be a scalar (1-dimensional vector).')
		if(as.integer(p.max)!=p.max)	stop('"p.max" must be integer.')
		if(p.max < 1)	stop('"p.max" must be a positive integer.')
		if(p.max >= 10000)	stop('"p.max" must be less than 10000.')

	  #check intercept and stand.coef
		candidate_TF <- c(TRUE, FALSE)
	 	if(sum(candidate_TF == intercept) != 1)	stop('intercept must be TRUE or FALSE.')
	 	if(sum(candidate_TF == stand.coef) != 1)	stop('stand.coef must be TRUE or FALSE.')




if(penalty=="enet" || penalty=="genet"){
		ex_para <- alpha
	}else{
		ex_para <- rep(0,2)
		ex_para[1] <- gamma
		ex_para[2] <- lambda
	}




	dfgps_result <- dfgps(X,y,penalty=penalty, ex_para=ex_para,  STEP=STEP, STEP.max=STEP.max, DFtype=DFtype, p.max= p.max)

	#	X0 <- scale(X) / sqrt(nrow(X)-1)
	#	y0 <- y-mean(y)
		
		
		
		#standardize
		#mean, varを計算
		
		#meanX <- apply(X,2,mean)
		#meanX_mat <- sweep(X, 2, meanX)
		#standardize_vec <- 1 / sqrt(apply(meanX_mat^2,2,sum))
		#standardize_mat <- matrix(rep(standardize_vec,nrow(X)),nrow(X),ncol(X),byrow=T)
		#meany <- mean(y)

		#X0 <- meanX_mat * standardize_mat
		#y0 <- y-meany


				if(missing(tau2)){
					Xstand = dfgps_result$Xstand
					ystand = dfgps_result$ystand
				if( nrow(Xstand) <= ncol(Xstand)){
						tau2 <- var(ystand)
					}else{
						tau2 <- sum((ystand-Xstand%*%dfgps_result$beta_OLS)^2) / (length(ystand) - ncol(Xstand))
					}
				}

	dfcp_result <- cp.dfgps(dfgps_result,tau2,intercept, stand.coef)
	dfaicc_result <- aicc.dfgps(dfgps_result,intercept, stand.coef)
	dfgcv_result <- gcv.dfgps(dfgps_result,intercept, stand.coef)
	dfbic_result <- bic.dfgps(dfgps_result,tau2,intercept, stand.coef)	
	ans <- list(stand.coef=stand.coef,intercept=intercept,dfgps_result=dfgps_result,dfcp_result=dfcp_result,dfaicc_result=dfaicc_result,dfgcv_result=dfgcv_result,dfbic_result=dfbic_result,penalty=penalty,alpha=alpha, gamma=gamma, lambda=lambda, tau2= tau2, call=match.call())
	class(ans) <- "msgps"
	ans
}

