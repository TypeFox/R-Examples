################################################################
##step 3a

cdlIBS=function(dat){

    if(have.cdlIBS(dat$snp.support)) {
      # cdlIBS info is already in snp.support. We only have to augment 
      # the cdlIBS probs by adding P(IBS=2|IBD=x) terms
      return(augment.cdlIBS(dat$snp.support))
    }
    # Otherwise, we need to calculate cdlIBS terms. 
    popsam=dat$subject.support$popsam
    snpobj0=dat$snp.data[popsam,]
    
    nbs=ncol(snpobj0)
    sumsnp <-  chopsticks::summary(snpobj0)
    
    nh <- 2*sumsnp$Calls # 2 times number of successful calls of each SNP
    X <- (nh/2)*(2*sumsnp$P.BB + 1*sumsnp$P.AB) # number of B alleles
    Y <- nh-X
    p <- X/nh
    q <- 1-p
    nbs=ncol(snpobj0)
    pi0z0=2*(p*q)^2*(X-1)*(Y-1)*nh^3/(X*Y*(nh-1)*(nh-2)*(nh-3))
    
    pi1z0=4*p^3*q*(X-1)*(X-2)*nh^3/(X*X*(nh-1)*(nh-2)*(nh-3))+
    4*p*q^3*(Y-1)*(Y-2)*nh^3/(Y*Y*(nh-1)*(nh-2)*(nh-3))
    
    pi2z0=p^4*(X-1)*(X-2)*(X-3)*nh^3/(X*X*X*(nh-1)*(nh-2)*(nh-3))+
    q^4*(Y-1)*(Y-2)*(Y-3)*nh^3/(Y*Y*Y*(nh-1)*(nh-2)*(nh-3))+
    4*(p*q)^2*(X-1)*(Y-1)*nh^3/(X*Y*(nh-1)*(nh-2)*(nh-3))
    pi0z1=rep(0,nbs)
    pi1z1=2*p^2*q*(X-1)*nh^2/(X*(nh-1)*(nh-2))+2*p*q^2*(Y-1)*nh^2/(Y*(nh-1)*(nh-2))
    pi2z1=p^3*(X-1)*(X-2)*nh^2/(X*X*(nh-1)*(nh-2))+
    q^3*(Y-1)*(Y-2)*nh^2/(Y*Y*(nh-1)*(nh-2))+
    p^2*q*(X-1)*nh^2/(X*(nh-1)*(nh-2))+p*q^2*(Y-1)*nh^2/(Y*(nh-1)*(nh-2))
    pi0z2=rep(0,nbs)
    pi1z2=rep(0,nbs)
    pi2z2=rep(1,nbs)
    
    outcdl=cbind(pi0z0,pi1z0,pi2z0,pi0z1,pi1z1,pi2z1,pi0z2,pi1z2,pi2z2)
    cnames<-c("pi0z0","pi1z0","pi2z0","pi0z1","pi1z1","pi2z1",
                  "pi0z2","pi1z2","pi2z2")
    colnames(outcdl)<-cnames
    
    return(outcdl)
}

have.cdlIBS<-function(ss) {
  if(!is.null(ss$pi0z0) && !is.null(ss$pi1z0) &&
     !is.null(ss$pi0z1) && !is.null(ss$pi1z1) &&
     !is.null(ss$pi0z2) && !is.null(ss$pi1z2) ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

add.cdlIBS<-function(ss,cdlibs) {
  # Add cdlibs to ss, or replace existing cdlibs columns in ss with 
  # those in cdlibs.
  if(!have.cdlIBS(ss)) {
    # add cdlibs: only need to save 2 columns out of each triplet
    return(cbind(ss,cdlibs[,-c(3,6,9)]))
  } else {
    # replace existing columns 
    ss$pi0z0<-cdlibs[,"pi0z0"]
    ss$pi1z0<-cdlibs[,"pi1z0"]
    ss$pi0z1<-cdlibs[,"pi0z1"]
    ss$pi1z1<-cdlibs[,"pi1z1"]
    ss$pi0z2<-cdlibs[,"pi0z2"]
    ss$pi1z2<-cdlibs[,"pi1z2"]
    return(ss)
  }
}

remove.cdlIBS<-function(ss) {
  ss$pi0z0<-ss$pi1z0<-ss$pi0z1<-ss$pi1z1<-ss$pi0z2<-ss$pi1z2<-NULL
  return(ss)
}

augment.cdlIBS<-function(ss) {
    # Have the cdlIBS info, just need to augment by P(IBS=2|IBD=x) terms
    out<-cbind(ss$pi0z0,ss$pi1z0,1-ss$pi0z0-ss$pi1z0,
               ss$pi0z1,ss$pi1z1,1-ss$pi0z1-ss$pi1z1,
               ss$pi0z2,ss$pi1z2,1-ss$pi0z2-ss$pi1z2)
    cnames<-cbind("pi0z0","pi1z0","pi2z0","pi0z1","pi1z1","pi2z1",
                  "pi0z2","pi1z2","pi2z2")
    colnames(out)<-cnames
    return(out)
}
