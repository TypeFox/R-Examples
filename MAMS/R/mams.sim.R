######################################################################
## 'mams.sim' evaluates 'generalsim' x times and computes average
######################################################################


mams.sim <- function(nsim=1000,nMat=matrix(c(44,88),nrow=2,ncol=5),u=c(3.068,2.169),l=c(0.000,2.169),pv=rep(0.5,4),ptest=1) { 

########################################################################################################
## 'generalsim' simulates the trial once. For general number of patients per arm per stage - 
## allocation given by the matrix R for active treatment arms. R[i,j]= allocation to stage i treatment j
## and vector r0 for control arm.
## treatment effect is specified by delta, delta0 and sig
##  Output is 
##  1) rejection yes/no
##  2) power yes/no
##  3) total sample size
#########################
  sim <- function(n,l,u,R,r0,delta,sig){
    J<-dim(R)[1]
    K<-dim(R)[2]
    Rdiff<-R-rbind(0,R[-J,])
    r0diff<-r0-c(0,r0[-J])

    ##############################################################################
    ## Create test statistics using independent normal increments in sample means:
    ##############################################################################

    mukhats<-apply(matrix(rnorm(J*K),J,K)*sig*sqrt(Rdiff*n)+Rdiff*n*matrix(delta,nrow=J,ncol=K,byrow=TRUE),2,cumsum)/(R*n)
    mu0hats<-cumsum(rnorm(J,0,sig*sqrt(r0diff*n)))/(r0*n)
    zks<-(mukhats-mu0hats)/(sig*sqrt((R+r0)/(R*r0*n)))

    ##############################################################################################################
    ## At each interim (j) determine whether or not treatment 1 has been declared significantly better than control 
    ## and also better than all other treatments:
    #############################################   

    eff<-0
    fut<-0
    ss<-0
    j<-1

    ### matrix of rejected hypotheses
    hmat <- matrix(0,nrow=J,ncol=K)
    remaining<-rep(T,K)
    while ((eff==0) && (fut==0) && (j<=J)) {
      ss<-sum((n*Rdiff[j,])[remaining])+n*r0diff[j]+ss
      eff<-(max(zks[j,remaining])>u[j])
      fut<-(max(zks[j,remaining])<l[j])
      if(any(zks[j,remaining]>u[j])){
        hmat[j,which((zks[j,]>u[j])&remaining)] <-1
      }       
      remaining<-(zks[j,]>l[j])&remaining
      j<-j+1

    }
    rej<-eff
    pow<-(remaining[1])&&eff
  
    if (pow){
      pow<-(which(zks[j-1,remaining]==max(zks[j-1,remaining]))==1)
    }

    if(all(remaining==0)){
      rem<-0
    }else{
      rem<-as.numeric(which(remaining))
    }

    return(list(rej=rej,pow=pow,ess=ss,hmat=hmat))
  }

  if(length(pv)!=(ncol(nMat)-1)) stop("Length of pv is not K")

  r0 <- nMat[,1]/nMat[1,1]
  R <-  nMat[,-1]/nMat[1,1]
  if(!is.matrix(R) && is.vector(R))  R <- t(as.matrix( nMat[, -1]/nMat[1, 1]))

  n <- nMat[1,1]
  deltas<-sqrt(2)*qnorm(pv)
  sig<-1

  reps<-sapply(rep(n,nsim),sim,l,u,R,r0,deltas,sig)

  ### power to reject any of the hypothesis corresponding to the treatments in ptest
  rej<-0
  for(i in 1:nsim){
    if(any(reps["hmat",i][[1]][,ptest]>0)){rej<-rej+1}
  }

  res <- NULL
  res$l <- l  
  res$u <- u
  res$n <- n
  res$K <- dim(R)[2]
  res$J <- dim(R)[1]

  res$rMat <- rbind(r0,t(R)) ## allocation ratios
  res$N <- sum(res$rMat[,res$J]*res$n) ## maximum total sample size


  res$nsim <- nsim
  res$typeI <- mean(unlist(reps["rej",]))
  res$power <- mean(unlist(reps["pow",]))
  res$ptest <- ptest
  res$prop.rej <- rej/nsim
  res$exss <- mean(unlist(reps["ess",]))

  class(res)<-"MAMS.sim"
  
  return(res)
}


#r<-mams.sim(10000,nMat=matrix(c(40,80),nrow=2,ncol=4),c(0,2.058),c(Inf,2.058),pv=c(0.65,0.55,0.55))
## type I error on treatment 1
#mean(r[2,])
## type I error overall 
#mean(r[1,])


