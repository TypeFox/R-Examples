laplace.fn <-
function(tune.vec,x.t,y.t,betahat.tm1,varbetahat.tm1) {

if (!is.matrix(x.t)) {
     dim(x.t) <- c(1,length(x.t))
	}
    Rhat.t <- varbetahat.tm1
    diag(Rhat.t) <- diag(Rhat.t) / tune.vec
	
    yhat.t <- dlogr.predict(x.t,betahat.tm1)
   
    Del1 <- t(x.t) %*% (y.t - yhat.t)
	Del2 <- -solve(Rhat.t) - (t(x.t) * matrix(rep(yhat.t*(1-yhat.t),dim(x.t)[2]),nrow=dim(x.t)[2],byrow=TRUE)) %*% x.t
	
    betahat.t <- betahat.tm1 - (solve(Del2) %*% Del1)
  #########################################################  
  #use rounding to deal with numerical issues with Rhat
  if(max(Rhat.t-t(Rhat.t))>0){p.theta <- dmnorm(t(betahat.t),t(betahat.tm1),round(Rhat.t,7))}
    else{p.theta <- dmnorm(t(betahat.t),t(betahat.tm1),Rhat.t)}
  #########################################################  

    

	p.y <- prod((exp(y.t*x.t %*%betahat.t))/(1 + exp(x.t %*% betahat.t)))
	
    return((((dim(x.t)[2])/2)*log(2*pi))+(.5*log(abs(det(ginv((1*Del2))))))+log(p.theta)+log(p.y))
}

