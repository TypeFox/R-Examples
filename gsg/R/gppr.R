gppr<-function(y,xterms,data,nterms=1,tol=0.001,
         gcvpen=1,maxit=50,family='binomial',
         max.terms=2){

    # family-specific functions
    fam<-get(family, mode = "function", envir = parent.frame())
    fam<-fam()
    g<-fam$linkfun; g.inv<-fam$linkinv; V<-fam$variance;          	
    g.prime<-function(x){1/V(x)}

    # initialize expectations
    response<-data[,as.character(y)]
    mu <- rep(mean(response), length(response))
    mu <- rep(0.5, length(response))
    deltaMu<-1; iter<-0;
    
    # the iterative re-weighting part
    gppr<-NULL
    f<-paste("Z","~",
         paste(xterms,collapse="+"),collapse="")
    while(deltaMu>=tol&iter<maxit){
        Z <- g(mu) + g.prime(mu) * (response - mu)
        cur.weights <- 1/(g.prime(mu)^2 * V(mu))
        data$Z<-Z
        gppr <- ppr(as.formula(f), weights=cur.weights,
                    sm.method="gcvspline",gcvpen=gcvpen,
                    nterms=nterms,data=data,
                    max.terms=max.terms,optlevel=3)
        eta <- predict(gppr,type='raw')

         muPrime<-mu; mu<-g.inv(eta);
         
         # deal with NaN that weights happen when mu%in%c(0,1)
         if(family=='binomial'){ 
         	boundary.tol<-10^(-5)
           mu[which(mu<boundary.tol)] <- boundary.tol
           mu[which(mu>(1-boundary.tol))] <- (1-boundary.tol)
         }
         deltaMu<-sum(abs(muPrime-mu))/sum(abs(muPrime))
         iter<-iter+1

    }
   
    output<-list(
      ppr=gppr
     ,family=fam
     ,iterations=iter
     ,data=data
     ,nterms=nterms
     ,tol=tol
     ,gcvpen=gcvpen
     ,maxit=maxit
     ,max.terms=max.terms
     ,y=y
     ,xterms=xterms
     ,formula=f
    )
    class(output)<-c("gppr")
    return(output)
}
