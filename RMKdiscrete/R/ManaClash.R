#None of these functions is important enough to have a C backend...
#TODO: handle p's equal to 0 for d functions, handle non-finite inputs better, 

#Bivariate PMF for damage done to opponent and you, either conditional or marginal w/r/t N:
dmanaclash.dmg <- function(x,y,N=NULL,pA=0.25,pB=0.25,pC=0.25,pD=0.25,log=FALSE){
  log <- log[1]
  if( any(pA<=0) | any(pB<=0) | any(pC<=0) | any(pD<=0) ){warning("NaNs produced")}
  if(!is.null(N)){
    #TODO: vectorize this, since it uses dbinom(), and apply() can be avoided
    inmtx <- cbind(x,y,N,pA,pB,pC,pD)
    return(apply(inmtx,1,function(row){
      .subdmg(x=row[1],y=row[2],N=row[3],pA=row[4],pB=row[5],pC=row[6],pD=row[7],log=log)})
    )
  }
  if(is.null(N)){
    inmtx <- cbind(x,y,pA,pB,pC,pD)
    return(apply(inmtx,1,function(row){
      .subdmg(x=row[1],y=row[2],N=NULL,pA=row[3],pB=row[4],pC=row[5],pD=row[6],log=log)})
    )
}}
.subdmg <- function(x,y,N,pA,pB,pC,pD,log){
  if(!is.finite(sum(x,y,N,pA,pB,pC,pD))){return(NA)}
  if(any(pA<=0,pB<=0,pC<=0,pD<=0)){return(NaN)}
  if( any(c(x,y)!=round(c(x,y))) ){return(ifelse(log==T,-Inf,0))}
  norm <- sum(pA,pB,pC,pD)
  pA <- pA/norm
  pB <- pB/norm
  pC <- pC/norm
  pD <- pD/norm
  r <- pA+pB+pC
  out <- 0
  if(!is.null(N)){ #CONDITIONAL on N
    out <- dbinom(x=x,size=N,prob=(pA+pC)/r,log=T) + dbinom(x=(x+y-N),size=x,prob=pA/(pA+pC),log=T) 
    return(ifelse(log==T,out,exp(out)))
  }
  if(is.null(N)){ #Marginal w/r/t N
    for(i in max(x,y):(x+y)){
      out <- out + 
        dmanaclash.xyN(x=x,y=y,N=i,pA=pA,pB=pB,pC=pC,pD=pD,log=F)
    }
    return(ifelse(log==T,log(out),out))
}}

#Trivariate joint PMF of damage to opponent, damage to self, and number of trials preceding double heads:
dmanaclash.xyN <- function(x,y,N,pA=0.25,pB=0.25,pC=0.25,pD=0.25,log=FALSE){
  log <- log[1]
  if( any(pA<=0) | any(pB<=0) | any(pC<=0) | any(pD<=0) ){warning("NaNs produced")}
  inmtx <- cbind(x,y,N,pA,pB,pC,pD)
  return(apply(inmtx,1,function(row){
    .subxyN(x=row[1],y=row[2],N=row[3],pA=row[4],pB=row[5],pC=row[6],pD=row[7],log=log)})
    )
}
.subxyN <- function(x,y,N,pA,pB,pC,pD,log){
  if(!is.finite(sum(x,y,N,pA,pB,pC,pD))){return(NA)}
  if(any(pA<=0,pB<=0,pC<=0,pD<=0)){return(NaN)}
  if( any(c(x,y,N)!=round(c(x,y,N))) ){return(ifelse(log==T,-Inf,0))}
  if(x>N | y>N){return(ifelse(log==T,-Inf,0))}
  if(any(c(x,y,N)<0)){return(ifelse(log==T,-Inf,0))}
  norm <- sum(pA,pB,pC,pD)
  pA <- pA/norm
  pB <- pB/norm
  pC <- pC/norm
  pD <- pD/norm
  out <- log(pD) - x*log(pB) + N*log(pB) - N*log(pA) + (x+y)*log(pA) - y*log(pC) + N*log(pC) + 
    lchoose(N,x) + lchoose(x,x+y-N)
  return(ifelse(log==T,out,exp(out)))
}


#PMF of net damage to opponent:
dmanaclash.net <- function(z,pA=0.25,pB=0.25,pC=0.25,pD=0.25,rel.eps=1e-8,log=FALSE){
  log <- log[1]
  if(any(length(z)>1,length(pA)>1,length(pB)>1,length(pC)>1,length(pD)>1)){
    if( any(pA<=0) | any(pB<=0) | any(pC<=0) | any(pD<=0) ){warning("NaNs produced")}
    inmtx <- cbind(z,pA,pB,pC,pD)
    return(apply(inmtx,1,function(row){
      dmanaclash.net(z=row[1],pA=row[2],pB=row[3],pC=row[4],pD=row[5],log=log)})
    )
  }
  else{
    if(!is.finite(sum(z,pA,pB,pC,pD,rel.eps))){return(NaN)}
    if(any(pA<=0,pB<=0,pC<=0,pD<=0)){return(NaN)}
    if(z!=round(z)){return(ifelse(log==T,-Inf,0))}
    norm <- sum(pA,pB,pC,pD)
    pA <- pA/norm
    pB <- pB/norm
    pC <- pC/norm
    pD <- pD/norm
    out <- 0
    if(sign(z)==1){x <- z; y <- 0}
    if(sign(z)==0){x <- y <- 0}
    if(sign(z)==-1){y <- -z; x <- 0}
    out <- dmanaclash.dmg(x=x,y=y,pA=pA,pB=pB,pC=pC,pD=pD,N=NULL,log=F)
    stopflag <- F
    while(stopflag==F){
      x <- x+1
      y <- y+1
      newterm <- dmanaclash.dmg(x=x,y=y,pA=pA,pB=pB,pC=pC,pD=pD,N=NULL,log=F)
      if(newterm/out <= rel.eps){stopflag <- T}
      out <- out+newterm
    }
    return(ifelse(log==T,log(out),out))
}}



rmanaclash <- function(n,pA=0.25,pB=0.25,pC=0.25,pD=0.25,N=NULL){
  n <- as.integer(n[1])
  out <- matrix(NA,nrow=n,ncol=3)
  colnames(out) <- c("x","y","N")
  pA <- rep(pA,length.out=n)
  pB <- rep(pB,length.out=n)
  pC <- rep(pC,length.out=n)
  pD <- rep(pD,length.out=n)
  norm <- pA+pB+pC+pD
  pA <- pA/norm
  pB <- pB/norm
  pC <- pC/norm
  pD <- pD/norm
  r <- pA+pB+pC
  if(!is.null(N)){
    N <- rep(N,length.out=n)
    inok <- rep(TRUE,length.out=n)
    if( any(pA<0) | any(pB<0) | any(pC<0) | any(pD<=0) | any(N!=round(N)) | any(!is.finite(pA+pB+pC+pD+N)) ){
      warning("NAs produced")
      inok[( pA<0 | pB<0 | pC<0 | pD<=0 | N!=round(N) | !is.finite(pA+pB+pC+pD+N) )] <- FALSE
      out[!inok,] <- NA
    }
    out[inok,1] <- rbinom(sum(inok),N[inok],(pA[inok]+pC[inok])/r[inok])
    out[inok,2] <- rbinom(sum(inok),out[inok,1],pA[inok]/(pA[inok]+pC[inok])) + 
      N[inok]-out[inok,1]
    out[inok,3] <- N[inok]
    }
  else{
    if( any(pA<0) | any(pB<0) | any(pC<0) | any(pD<=0) ){warning("NAs produced")}
    for(i in 1:n){
      if(pA[i]<0 | pB[i]<0 | pC[i]<0 | pD[i]<0){out[i,] <- NA; next}
      #TOmaybeDO: probably a faster way to generate random draws can be used, using rgeom() and rbinom().
      x <- 0
      y <- 0
      N <- 0
      stoptoss <- FALSE
      while(stoptoss==F){
        U <- runif(1,0,1)
        if(U<=(pA[i]+pB[i])){
          if(U<=pA[i]){x <- x+1; y <- y+1; N <- N+1; next}
          else{y <- y+1; N <- N+1; next}
        }
        else{
          if(U<=(pA[i]+pB[i]+pC[i])){x <- x+1; N <- N+1; next}
          else{stoptoss <- TRUE}
        }  
#         if(U<=pA[i]){x <- x+1; y <- y+1; N <- N+1; next}
#         if(U>pA[i] & U<=(pA[i]+pB[i])){y <- y+1; N <- N+1; next}
#         if(U>(pA[i]+pB[i]) & U<=(pA[i]+pB[i]+pC[i])){x <- x+1; N <- N+1; next}
#         if(U>(pA[i]+pB[i]+pC[i])){stoptoss <- TRUE}
      }
      out[i,] <- c(x,y,N)
    }
  }
  return(out)
}
