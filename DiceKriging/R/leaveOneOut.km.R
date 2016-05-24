`leaveOneOut.km` <-
function(model, type, trend.reestim=FALSE) {
	
  if (length(model@noise.var)>0) {
    stop("At this stage, leave-one-out is not implemented for noisy observations")
  }

  case1 <- ((type=="SK") && (!trend.reestim))
  case2 <- ((type=="UK") && (trend.reestim))
  analytic <- (case1 | case2)
  # If analytic, we can use Dubrule's formulae

  X <- model@X
	y <- model@y
	T <- model@T
	n <- nrow(as.matrix(X))
	yhat <- sigma2 <- matrix(NA,n,1) 

	beta <- model@trend.coef
	F <- model.matrix(model@trend.formula, data=data.frame(X,y))

	if (!analytic)	{	
	  
    C <- t(T)%*%T            # should be improved in future versions  
	  
		for (i in 1:n) {
      
			F.i <- F[i,]
			y.but.i <- y[-i,]
			F.but.i <- F[-i,]
			C.but.i <- C[-i,-i]
			c.but.i <- C[-i,i]
			T.but.i <- chol(C.but.i)					
			x <- backsolve(t(T.but.i), y.but.i, upper.tri=FALSE)
			M <- backsolve(t(T.but.i), F.but.i, upper.tri=FALSE)
			M <- as.matrix(M)
      
      if (trend.reestim) {    # reestimation of beta
        l <- lm(x ~ M - 1)
        beta <- as.matrix(l$coef, ncol=1)
      }
			z <- x - M%*%beta
			Tinv.c <- backsolve(t(T.but.i), c.but.i, upper.tri=FALSE)      # only a vector in this case
			y.predict.complement <- t(Tinv.c) %*% z
			y.predict.trend <- F.i %*% beta
			
			y.predict <- y.predict.trend + y.predict.complement
			yhat[i] <- y.predict
			
			sigma2.1 <- crossprod(Tinv.c)

      total.sd2 <- C[i,i]
      sigma2[i] <- total.sd2 - sigma2.1
			
      if (type=="UK"){
        T.M <- chol(t(M)%*%M)
			  sigma2.mat <- backsolve(t(T.M), t(F.i - t(Tinv.c)%*%M) , upper.tri=FALSE)
			  sigma2.2 <- apply(sigma2.mat, 2, crossprod)
			  sigma2[i] <- sigma2[i] + sigma2.2
			}
      
      sigma2[i] <- max(sigma2[i], 0)
			sigma2[i] <- as.numeric(sigma2[i])
			
		}	# end 'for' loop
	
  }	else {   # fast computation
    Cinv <- chol2inv(T)       # cost : n*n*n/3
    if (trend.reestim & (type=="UK")) {
      M <- model@M              # recall : M = inv(t(T))*F  (cost to recompute : n*n*p)
      Cinv.F <- Cinv %*% F      # cost : 2*n*n*p
      T.M <- chol(crossprod(M))   # cost : p*p*p/3, neglected
      aux <- backsolve(t(T.M), t(Cinv.F), upper.tri=FALSE)   # cost : p*p*n, neglected
      Q <- Cinv - crossprod(aux)    # cost : 2*n*n*(p-1/2)
      Q.y <- Q%*%y
      # Remark:   Q <- Cinv - Cinv.F %*% solve(t(M)%*%M) %*% t(Cinv.F)   # direct (not so bad actually)
    } else if ((!trend.reestim) & (type=="SK")){
      Q <- Cinv
      Q.y <- Q%*%(y-F%*%beta)
    } else {      
      stop("This case is not implemented yet")
    }
    sigma2 <- 1/diag(Q)
    epsilon <- sigma2 * (Q.y)  # cost : n, neglected 
    yhat <- as.vector(y - epsilon)
  }
  
  
  return(list(mean=as.numeric(yhat), sd=as.numeric(sqrt(sigma2))))
}