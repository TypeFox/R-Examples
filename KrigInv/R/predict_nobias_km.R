
predict_nobias_km <- function(object, newdata, type="UK", se.compute=TRUE, cov.compute=FALSE,low.memory=FALSE,...) {
	#newdata : n x d

	X <- object@X
	y <- object@y
	
	newdata <- as.matrix(newdata)
	dim.newdata <- dim(newdata)
	m <- dim.newdata[1]
	d.newdata <- dim.newdata[2]
	
	if (!identical(d.newdata, object@d)) stop(paste("newdata must have the same numbers of columns (",d.newdata,") than the experimental design (",object@d,")"))
	if (!identical(colnames(newdata), colnames(X))) colnames(newdata) <- colnames(X)
	
	T <- object@T
	z <- object@z
	M <- object@M
	beta <- object@trend.coef
	F.newdata <- model.matrix(object@trend.formula, data=data.frame(newdata))
	y.predict.trend <- F.newdata%*%beta
	
	c.newdata <- covMat1Mat2(X1=X, X2=newdata, object=object@covariance, nugget.flag=object@covariance@nugget.flag)    # compute c(x) for x = newdata ; remark that for prediction (or filtering), cov(Yi, Yj)=0  even if Yi and Yj are the outputs related to the equal points xi and xj.
	
	Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri=FALSE)
	y.predict.complement <- crossprod(Tinv.c.newdata,z)
	y.predict <- y.predict.trend + y.predict.complement
	y.predict <- as.numeric(y.predict)
	
  if(low.memory){ 
    output.list <- list(mean = y.predict)
  } else {
    output.list <- list(mean = y.predict, c=c.newdata, Tinv.c=Tinv.c.newdata,F.newdata=F.newdata)
  }
  
	if (se.compute==TRUE) {		
		
		s2.predict.1 <- apply(Tinv.c.newdata, 2, crossprod)         # compute c(x)'*C^(-1)*c(x)   for x = newdata
			
		if (object@covariance@nugget.flag) {
			total.sd2 <- object@covariance@sd2 + object@covariance@nugget
		} else total.sd2 <- object@covariance@sd2
		
		
		if (type=="SK") {
			s2.predict <- pmax(total.sd2 - s2.predict.1, 0)
			s2.predict <- as.numeric(s2.predict)
			q95 <- qnorm(0.975)
		}
		else if (type=="UK") {
			T.M <- chol(crossprod(M))   # equivalently : qrR <- qr.R(qr(M))
			s2.predict.mat <- backsolve(t(T.M), t(F.newdata - crossprod(Tinv.c.newdata,M)) , upper.tri=FALSE)
			
			s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
			s2.predict <- pmax(total.sd2 - s2.predict.1 + s2.predict.2, 0)
			s2.predict <- as.numeric(s2.predict)


			q95 <-  qnorm(0.975)
			
		}
		
		lower95 <- y.predict - q95*sqrt(s2.predict)
		upper95 <- y.predict + q95*sqrt(s2.predict)
		output.list$sd <- sqrt(s2.predict)
		if(!low.memory) output.list$lower95 <- lower95
		if(!low.memory) output.list$upper95 <- upper95
	}
	
	if (cov.compute==TRUE) {		
		
		if (object@covariance@nugget.flag) {
			total.sd2 <- object@covariance@sd2 + object@covariance@nugget
		} else total.sd2 <- object@covariance@sd2
		
		C.newdata <- covMatrix(X=newdata, object=object@covariance)[[1]]
		cond.cov <- C.newdata - crossprod(Tinv.c.newdata)
		
		if (type=="UK") {	
			T.M <- chol(t(M)%*%M)
			s2.predict.mat <- backsolve(t(T.M), t(F.newdata - t(Tinv.c.newdata)%*%M), upper.tri=FALSE)
			cond.cov <- cond.cov + crossprod(s2.predict.mat)
		}
		
		output.list$cov <- cond.cov
		
	}
	return(output.list)
}
