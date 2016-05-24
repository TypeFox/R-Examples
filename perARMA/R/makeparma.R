makeparma <-
function(n,phi,theta,del,nprep) {

    Tp=nrow(phi)
    p=ncol(phi)
    Tq=nrow(theta)
    q=ncol(theta)
    Td=nrow(del)
    dum=ncol(del)

     if ((p & (Tp!=Td)) | (q & (Tq!=Td)))  {stop("period incompatibility")}
    T=Td

    a=cbind(matrix(1,T,1),-phi)    
    b=cbind(del,theta)               
      if (nargs()<5) {nprep=50}
      if (p==0) {nprep=0}

       ntot=n+nprep*T                              
      xi=rnorm(ntot)
      xi=as.matrix(xi)
                                        
     parmafil_out<-parmafil(b,a,xi)
     y=parmafil_out$y
     y=y[(nprep*T+1):ntot]                     

      result = list(y=y)   
      class(result) = "makeparma"
      result
}

