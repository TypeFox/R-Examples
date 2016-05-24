dbglm.iteration <- function(y, mu, weights,  nobs, eta, Delta, method, offset, n, eff.rank = NULL, rel.gvar, dev.resids,
                                         aic, mu.eta, valideta, validmu, family,variance, linkinv, problem.links,
                                         eps1,eps2,maxiter)
{
  # Deviance Dev
   Dev      <-sum(dev.resids(y,mu,weights))
   dblm_aux <- NULL
   # Step i, i=iter, 1<=iter<=maxiter
   iter <- 0
   while ( TRUE ){
     iter <- iter + 1
     mu0 <- mu
     Dev0 <- Dev
     dblm0 <- dblm_aux
     eta0 <-eta

     # gdmu and wdmu, using the derivative function of eta
     gdmu<- 1/mu.eta(eta)
     wdmu<- mu.eta(eta)

     # the weights must be > 0 and wdmu must be != 0 to fit the dblm.
     good <- weights > 0
     if (any(is.na(wdmu[good])))
       stop("NAs in d(mu)/d(eta)")
     good <- (weights > 0) & (wdmu != 0)
     # if all the observations are not good -->  break the program
     if (all(!good)) {
        stop("no observations informative at iteration ",iter)
     }
     # varmu = variance only with good observations.
     varmu <- variance(mu)[good]
     if (any(is.na(varmu)))
       stop("NAs in V(mu)")
     if (any(varmu == 0))
       stop("0s in V(mu)")

     # new y, to fit the dblm. z only has the good observations
     z <- (eta-offset)[good] + (y-mu)[good]/wdmu[good]

     # new_weights, see pag 43 of generalized linear models.
     new_weights<- (weights[good] * wdmu[good]^2)/varmu
     Delta_aux<-Delta[good,good]

     # make the dblm only with the valid observations (weights > 0)
     class(Delta_aux)="D2"

     if (!is.null(eff.rank))
      method="eff.rank"
     else
      method="rel.gvar"

    if (!is.null(eff.rank))
       eff.rank_aux<-min(length(z)-1,eff.rank)
    else
       eff.rank_aux<- NULL

     if (length(z)<2)
      stop(paste("there are only one observation in the dblm iteration ",iter,". The weights are negative. Try to use other link"))


     dblm_aux <- dblm.D2(D2=Delta_aux,y=z,weights=new_weights,method=method,
                    rel.gvar=rel.gvar,eff.rank=eff.rank_aux)

     Hhat <- dblm_aux$H

     # eta[good] = fitted values of the dblm (only with good observations)
     eta[good]<-dblm_aux$fitted.values

     # predict the values of the observations with negative weights or a wdmu=0
     if (any(!good)){
       Hhat=matrix(0,ncol=nobs,nrow=nobs)

       Delta_pr<-Delta[!good,good]
       if (is.null(nrow(Delta_pr)))
        Delta_pr <- t(as.matrix(Delta_pr))

       class(Delta_pr)<-"D2"
       eta[which(!good)]<-predict(dblm_aux,Delta_pr,type="D2")
       Hhat[good,good]<-dblm_aux$H
     }

     # link: 1. logit, 2. logarithmic, 3. Identity, 4. mu-1, 5. mu-2     
     mu<-linkinv(eta <- eta+ offset) # Aquí he canviat
     Dev <- sum(dev.resids(y,mu,weights))

     if (is.finite(any(mu<0)))
      neg.mu<-any(mu<0)
     else
      neg.mu<-TRUE

     if (!is.finite(Dev) && neg.mu && any(problem.links==family$link)){
      stop(gettextf("The linear predictor is negative and some components of 'mu' are < 0. The Deviance in the iteration %i is NULL,
                      Try to change the link to any one different than: 'inverse', 'identity', 'logit' or '1/mu^2'",iter))
     }

     if (!is.finite(Dev))
      stop("The deviance is Null. Try to change the link function.")

     # check if the the relative error of mu and Dev is small enough
     rel_incr_mu <- sum(abs(mu-mu0))/sum(abs(mu))
     rel_incr_Dev <- abs(Dev - Dev0)/(0.1 + abs(Dev))

     if (rel_incr_Dev < eps1){
       convcrit <- "DevStat"
       break
     }
     if (rel_incr_mu < eps2){
       convcrit <- "muStat"
       break
     }
     if (iter>=maxiter){
       convcrit <- "MaxIter"
       break
     }
   }

   # residuals, residuals degree of freedom, null deviance and aic
   wtdmu         <-sum(weights*y)/sum(weights)
  
  null.deviance <- eval(call("glm.fit", 
     x = rep(1,length(y)), y = y, weights = weights, 
     offset = offset, family = family, intercept = TRUE))$deviance
     
  # null.deviance <-sum(dev.resids(y,wtdmu,weights))
   n.ok          <-nobs-sum(weights==0)
   df.null       <-n.ok-1
   df.residual   <-df.null-dblm_aux$eff.rank
   aic.model     <-aic(y,nobs,mu,weights,Dev)+2*(dblm_aux$eff.rank+1)
   bic.model     <- aic(y, nobs, mu, weights, Dev) + log(nobs) * (dblm_aux$eff.rank+1)
   gcv.model     <- n.ok^2* (Dev/sum(weights)) / (n.ok - sum(diag(Hhat)))^2
   residuals     <- (y-mu)/mu.eta(eta)

   
   return(list(aic.model = aic.model, bic.model = bic.model, gcv.model = gcv.model,
               df.residual = df.residual, residuals=residuals, fitted.values=mu,
               family=family, deviance=Dev, null.deviance=null.deviance,iter=iter,
               weights=new_weights,prior.weights=weights,
               varmu=variance(mu), dev.resids=dev.resids(y,mu,weights), 
               df.null=df.null,y=y,convcrit=convcrit,H=Hhat,
               rel.gvar=dblm_aux$rel.gvar,eff.rank=dblm_aux$eff.rank, dblm_aux=dblm_aux))
}