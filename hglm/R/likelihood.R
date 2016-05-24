`likelihood` <- function(hglm.obj, y, X, Z, weights = NULL,family=NULL) {
	if(!is.null(hglm.obj$likelihod)) return(hglm.obj$likelihod)
  hgcall<-hglm.obj$call
	familyM<-hgcall$family
	if(is.null(familyM)){
  if(!is.null(family)){
      familyM<-family
      if (is.character(familyM)) familyM <- get(familyM)
      if (is.function(familyM)) familyM <- eval(familyM)
        }else{
      familyM<-gaussian() }
  }else{
   familyM <- eval(familyM)
  }
  if(((familyM$family[1]=="binomial") | (familyM$family[1]=="poisson"))& (is.null(hgcall$fix.disp))) stop("Likelihoods are not available for quasi models")
  n <- length(y)
	if (is.null(hglm.obj$disp.fv)) hglm.obj$disp.fv <- rep(hglm.obj$varFix, n)
	if (is.null(weights)) weights <- rep(1,n)
	 glm.phi <- as.numeric(hglm.obj$disp.fv)/weights
	const <- switch(familyM$family[1],
			gaussian = -0.5*sum(log(2 * pi*glm.phi)),
			poisson =  -1*sum(y - y * ifelse(y > 0,log(y), 0) + lgamma(y + 1)),
			Gamma = -1*sum(log(y)+1/glm.phi+lgamma(1/glm.phi)+log(glm.phi)/glm.phi),
			inverse.gaussian = -0.5* sum(log(2 * pi * y^3)),
			binomial = 0
	)
	
	#linear.pred <- family$linkfun(hglm.obj$fv)
	#linear.pred_disp <- hglm.obj$disp.fv
	#h.dglm <- dglm(y ~ 0 + offset(linear.pred), ~ 0 + offset(linear.pred_disp), family = family, dlink="identity")
	#cond.lik <- summary(h.dglm)$m2loglik
   logdet <- function(M){
   eigA <- eigen(M, only.values=TRUE)$values
  condA <- min(eigA)/max(eigA)
  if (condA<1e-8) warning("The Hessian used for computing pbvh is ill-conditioned.")
  return(sum(log(eigA)))
  }
	eta.i<- familyM$linkfun(hglm.obj$fv[1:n])
	dmu_deta <- familyM$mu.eta(eta.i)
	d <- hglm.obj$dev[1:n]
	cond.lik <- const -0.5* sum(d/glm.phi)
	if((familyM$family[1]=="binomial") & (!is.logical(all.equal(weights, rep(1,n))))) cond.lik <- sum(dbinom(x=round(y*weights,0),size=weights,prob=hglm.obj$fv[1:n],log=TRUE))
	p <- sum(hglm.obj$hv[1:n])
	cAIC <- -2*cond.lik + 2*p
	v <- hglm.obj$ranef
	RandC<-hglm.obj$RandC
	nRand<-length(RandC)
	RanF<-hglm.obj$call.rand.family
	if(!is.null(RanF)) RanF<-eval(RanF)
	hlik<-cond.lik
	vs<-1
	vend<-0
  for(i in 1:nRand){
	if(is.null(RanF[[i]])){familyR<-gaussian()
  }else{
    if(class(RanF)=="family"){
     familyR<-RanF
    }else{
  familyR<-eval(RanF[[i]]) }
  }
  vend<-vend+RandC[i]
  u.i<-v[vs:vend]
  v.i<-familyR$linkfun(u.i)
	du_dv<- familyR$mu.eta(v.i)
	RanDisp<-hglm.obj$phi[vs:vend]
	dh2_dv2<-(du_dv^2/familyR$variance(u.i))*(1/RanDisp)
	lfv<-switch(familyR$family[1],
   gaussian = sum(dnorm(u.i,0,sd=sqrt(RanDisp),log=TRUE)),
		Beta =  sum(dbeta(x=u.i,shape1=1/(2*RanDisp),shape2=1/(2*RanDisp),log=TRUE)) ,
		Gamma = sum(dgamma(x=u.i,shape=1/RanDisp,scale=RanDisp,log=TRUE)),
		GAMMA = sum(dgamma(x=u.i,shape=1/RanDisp,scale=RanDisp,log=TRUE)),
		CAR = sum(dnorm(x=(t(familyR$Dvec)%*%u.i),0,sd=sqrt(hglm.obj$CAR.tau/(1-hglm.obj$CAR.rho*familyR$Deigen)),log=TRUE)),
		SAR = sum(dnorm(x=(t(familyR$Dvec)%*%u.i),0,sd=sqrt(hglm.obj$SAR.tau/(1-hglm.obj$SAR.rho*familyR$Deigen)^2),log=TRUE))
		#Inverse-Gamma is not implemented ######	
  )
  if(i==1){
  w2<-dh2_dv2
  }else{
  w2<-c(w2,dh2_dv2)
  }
	hlik <- hlik +lfv +sum(log(abs(du_dv)))
	vs<-vs+RandC[i]
  }
  W1 <- diag(as.numeric(dmu_deta^2/(glm.phi*familyM$variance(hglm.obj$fv[1:n]))))
  W2 <- diag(w2)
	H <- t(Z)%*%W1%*%Z+W2
	A <- rbind(cbind(t(X)%*%W1%*%X,t(X)%*%W1%*%Z),cbind(t(Z)%*%W1%*%X,H))
	pvh <- hlik -0.5* logdet((H/(2*pi)))# hlik  + log(abs(det(H/(2*pi))))
	## xia 2014-01-21: pvh calculations bypass solve(H) and solve(W2) to avoid numerical problems while hglm.obj$varRanef is very small
	pbvh <- hlik -0.5* logdet(A/(2*pi))
	list(hlik = as.numeric(hlik), pvh = as.numeric(pvh), pbvh = as.numeric(pbvh), cAIC = cAIC)
}
