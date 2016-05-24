fnnls <- 
  function(XtX,Xty,ntol=NULL){
    # Fast Non-Negative Least Squares algorithm based on 
    #   Bro, R., & de Jong, S. (1997) A fast non-negativity-constrained 
    #   least squares algorithm. Journal of Chemometrics, 11, 393-401.
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    ### initialize variables
    pts <- length(Xty)
    if(is.null(ntol)){
      ntol <- 10*(.Machine$double.eps)*max(colSums(abs(XtX)))*pts
    }
    pvec <- matrix(0,1,pts)
    Zvec <- matrix(1:pts,pts,1)
    beta <- zvec <- t(pvec)
    zz <- Zvec
    wvec <- Xty - XtX%*%beta
    
    ### iterative procedure
    iter <- 0    
    itmax <- 30*pts
    
    # outer loop
    while(any(Zvec>0) && any(wvec[zz]>ntol)) {
      
      tt <- zz[which.max(wvec[zz])]
      pvec[1,tt] <- tt
      Zvec[tt] <- 0
      pp <- which(pvec>0)
      zz <- which(Zvec>0)
      nzz <- length(zz)
      zvec[pp] <- smpower(XtX[pp,pp],-1)%*%Xty[pp]
      zvec[zz] <- matrix(0,nzz,1)
      
      # inner loop
      while(any(zvec[pp]<=ntol) &&  iter<itmax) {
        
        iter <- iter + 1
        qq <- which((zvec<=ntol) & t(pvec>0))
        alpha <- min(beta[qq]/(beta[qq]-zvec[qq]))
        beta <- beta + alpha*(zvec-beta)
        indx <- which((abs(beta)<ntol) & t(pvec!=0))
        Zvec[indx] <- t(indx)
        pvec[indx] <- matrix(0,1,length(indx))
        pp <- which(pvec>0)
        zz <- which(Zvec>0)
        nzz <- length(zz)
        if(length(pp)>0){
          zvec[pp] <- smpower(XtX[pp,pp],-1)%*%Xty[pp]
        }
        zvec[zz] <- matrix(0,nzz,1)      
        
      } # end inner loop
      
      beta <- zvec
      wvec <- Xty - XtX%*%beta
      
    } # end outer loop
    
    beta
    
  }