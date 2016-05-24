#########################################
## Modelos  Mistos Censurados - Normal ##
#########################################

## Momentos da Normal Truncada ##

MomemNT<-function(u=c(0,0),S=diag(2),qc=c(1,2)) {

  nic=length(u)

  if (nic==1) {

    qq <- (1/sqrt(S))*(-qc+u)
    R<-1
    alpha <- pnorm(-qq)

    dd <- dnorm(-qq)
    H <- qq*dd
    EX <- (1/alpha)*dd   # a vector with a length of nic
    EXX <- 1+1/alpha*H
    varX <- EXX-EX^2
    Eycens <- -sqrt(S)*EX+u
    varyic<- varX*S
    E2yy<-varyic+Eycens^2

  }

  else {

    qq <- diag(1/sqrt(diag(S)))%*%(-qc+u)
    R <-  diag(1/sqrt(diag(S)))%*%S%*%diag(1/sqrt(diag(S)))
    alpha <- pmvnorm(upper=as.vector(-qq), corr=R)
    dd <- rep(0, nic)   #derivative vector

    for (j in 1:nic){
      V <- R[-j, -j, drop=F]-R[-j,j, drop=F]%*%R[j,-j, drop=F]
      nu <- -qq[-j]+R[-j,j, drop=F]%*%qq[j]
      dd[j] <- dnorm(-qq[j])*pmvnorm(upper=as.vector(nu), sigma=V)
    }

    H <- matrix(rep(0, nic*nic), nrow=nic)
    RH <- matrix(rep(0, nic*nic), nrow=nic)

    if(nic==2)     {
      H[1,2] <- H[2,1] <- dmvnorm(-qq[c(1, 2)],sigma=matrix(c(1, R[1,2], R[2,1], 1), nrow=2))
      #sigma==R since qq is standardized
      RH[1,2] <- RH[2,1] <- R[1,2]*H[1,2]
    }

    else {
      for( s in 1:(nic-1)){
       for (t in (s+1):nic){
        invR <- solve(R[c(s,t), c(s,t), drop=F])
        nu <- -qq[-c(s,t)]+R[-c(s,t), c(s,t), drop=F]%*%invR%*%qq[c(s,t),,drop=F]
        V <-  R[-c(s,t), -c(s,t), drop=F]- R[-c(s,t), c(s,t), drop=F]%*%invR%*%R[c(s,t), -c(s,t), drop=F]
        H[s,t] <- H[t,s] <- pmvnorm(upper=as.vector(nu), sigma=V)*dmvnorm(-qq[c(s, t)],sigma=matrix(c(1, R[s,t], R[t,s], 1), nrow=2))
        RH[s,t] <- RH[t,s] <- R[s,t]*H[s,t]
       }
      }
    }

    h <- qq*dd-apply(RH, 1, sum)
    diag(H) <- h
    EX <- (1/alpha)*R%*%dd   # a vector with a length of nic
    EXX <- R+1/alpha*R%*%H%*%R
    varX <- EXX-EX%*%t(EX)
    Eycens <- -diag(sqrt(diag(S)))%*%EX+u
    varyic <- diag(sqrt(diag(S)))%*%varX%*%diag(sqrt(diag(S)))
    E2yy <- varyic+Eycens%*%t(Eycens)

  }

  return(list(Ey=Eycens,Eyy=E2yy,Vary=varyic))

}
