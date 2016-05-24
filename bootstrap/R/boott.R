"boott"<- function(x,theta,...,sdfun=sdfunboot,nbootsd=25,nboott=200,
		   VS=FALSE, v.nbootg=100,v.nbootsd=25,v.nboott=200,
		   perc=c(.001,.01,.025,.05,.10,.50,.90,.95,.975,.99,.999)){

    call <- match.call()
    sdfunboot <- function(x,nboot,theta,...){
        n <- length(x)
        junk <- matrix(sample(x,size=n*nboot,replace=TRUE),nrow=nboot)
        return(sqrt(var(apply(junk,1,theta,...))))
    }

    thetahat<- theta(x,...)
    n <- length(x)
    if(!VS) {sd <- sdfun(x,nbootsd,theta,...)} else {sd <- 1}
    
    if(VS){
        xstar <- matrix(sample(x,size=n*v.nbootg,replace=TRUE),
                        nrow=v.nbootg)
        thetastar0 <- apply(xstar,1,theta,...)
        sdstar0 <- apply(xstar,1,sdfun,v.nbootsd,theta,...)
        o <- order(thetastar0)
        thetastar0 <- thetastar0[o]
        sdstar0 <- sdstar0[o]
        
        temp <- lowess(thetastar0,log(sdstar0))$y
        
        sdstar0 <- exp(temp)
        invsdstar0 <- 1/sdstar0
        g <- ctsub(thetastar0,invsdstar0,thetastar0)
        g <- (g-mean(g))/sqrt(var(g))
        g <- g*sqrt(var(thetastar0))+mean(thetastar0)
    }

    if(!VS) { thetastar0 <- NULL; g  <-  NULL}

    if(!VS) {
        xstar <- matrix(sample(x,n*nboott,replace=TRUE),nrow=nboott)
    } else {
        xstar <- matrix(sample(x,n*v.nboott,replace=TRUE),nrow=v.nboott)
    }
    thetastar <- apply(xstar,1,theta,...)
    gthetastar <- rep(0,length(thetastar))
    
    if(VS) {
        gthetahat <- yinter(thetastar0,g,thetahat)
    } else {
        gthetahat <- thetahat
    }
    
    if(VS){
        for(i in 1:length(thetastar)){
            gthetastar[i] <- yinter(thetastar0,g,thetastar[i])
        }
    }
    else {
        gthetastar <- thetastar
    }
  
    if(!VS) {
        sdstar <- apply(xstar,1,sdfun,nbootsd,theta,...)
    } else {
        sdstar <- 1
    }
  
    tstar <- sort( (gthetastar-gthetahat)/sdstar)[length(gthetastar):1]

    ans <-  gthetahat-sd*tstar

    if(VS){
        for(i in 1:length(ans)) {
            ans[i] <- xinter(thetastar0,g,ans[i])
        }
    }

    o <- trunc(length(ans)* perc)+1
  
    ans1 <- matrix(ans[o],nrow=1)

    dimnames(ans1) <- list(NULL,perc)

    return(list(confpoints=ans1,
                theta=thetastar0,
                g=g,
                call=call))
}
