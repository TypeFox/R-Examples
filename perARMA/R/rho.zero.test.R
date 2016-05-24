rho.zero.test <-
function(rhovec,nvec) {

     nrho=length(rhovec)
     nn=length(nvec)

    if (nrho != nn) {stop(" nrho != nn ")}
    ind1=which(abs(rhovec-1) < 2*.Machine$double.eps)
     rhovec[ind1]=1
     zval<-matrix(1,1,nrho)

     good=which(rhovec!=1)
     ngood=length(good)
      zval[good]=0.5*log((1+rhovec[good])/(1-rhovec[good]))   
      zvar=1/(nvec[good]-3)   
      zstd=sqrt(zvar)         
      zzval=zval/zstd         
      mzscore=mean(zzval)*sqrt(ngood)
   if (mzscore)
   { pv=2*(1-pnorm(abs(mzscore),0,1))  
    } else {
      pv=1}

    result = list(pv=pv,mzscore=mzscore,ngood=ngood)
    class(result) = "rho.zero.test"
    result
}

