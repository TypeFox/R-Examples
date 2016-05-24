rho.constant.test <-
function(rhovec,nvec){

       nrho=length(rhovec)
       nn=length(nvec)

       if (nrho != nn)
        {stop("nrho != nn")}

       zval<-matrix(1,1,nrho)
       good=which(rhovec!=1)
       zval[good]=0.5*log((1+rhovec[good])/(1-rhovec[good]))   

       zetahat=mean(zval) 
       zvar=1/(nvec[good]-3)
       zval=zval-zetahat
       chi2=base::sum((zval*zval)/zvar)   
       ngood=length(good)

      if (chi2){
         pv=1-pchisq(chi2,ngood-1)
          } else {pv=1}

    result = list(pv=pv,chi2=chi2,ngood=ngood)
    class(result) = "rho.constant.test"
    result
}

