f.NI <- function(coords,params,rel.tol,abs.tol,K,cov.f,info)
  {
    a <- info$rowwidth
    b <- info$colwidth
    ax <- coords[1]
    bx <- coords[2]
    i <- coords[3]
    j <- coords[4]
    ## 5 
    ## 6
    ## 7 - which rows to evaluate - evalFactor
    ## 8 - lengths
    ## 9 - indicator for which rows to evaluate analytic results
    minD <- coords[10]
    maxD <- coords[11]

    if(minD < maxD) {
      int.res <- try(integrate(cov.f,lower=minD, upper=maxD,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j,rel.tol=rel.tol,abs.tol=abs.tol,params=params,K=K,stop.on.error=F))
    } else { int.res <- list(value=0) }
    
    if(class(int.res)=="try-error") {

      cat("integration fail for plot",i,j,"with limits of integration",minD,maxD,"corresponding to",integrate(f,lower=minD, upper=maxD,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j,rel.tol=rel.tol,abs.tol=abs.tol)$value,"of the density and over this region the covariance function values are between",cov.f(minD,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j,params=params,K=K),"and",cov.f(maxD,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j,params=params,K=K),"\n")

      if(FALSE) {
        par(mfrow=c(3,2))
        gVals <- seq(0,1,length=101)
        gVals <- c(gVals^2,1-gVals^2)
        gVals <- gVals*(maxD-minD)+minD
        gVals <- sort(gVals)
        plot(gVals,K(gVals,params=params),xlab="",ylab="",main="Cov Function K")
        plot(gVals,f(gVals,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j),xlab="",ylab="",main="Density f")
        plot(gVals,cov.f(gVals,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j,params=params,K=K),xlab="",ylab="",main="K times f")
        gVals <- seq(0,1,length=101)
        gVals <- c(gVals^2,1-gVals^2)
        gVals <- gVals*(maxD-minD+4)+minD-4
        gVals <- sort(gVals)
        plot(gVals,K(gVals,params=params),xlab="",ylab="",main="Cov Function K, larger range")
        plot(gVals,f(gVals,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j),xlab="",ylab="",main="Density f, larger range")
        plot(gVals,cov.f(gVals,rw=a,cw=b,ax=ax,bx=bx,i=i,j=j,params=params,K=K),xlab="",ylab="",main="K times f, larger range")
        locator(1)
      }
      
      int.res <- list(NULL)
      int.res$value <- 0
    }
    int.res$value
  }
