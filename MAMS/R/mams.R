mams <- function(K=4, J=2, alpha=0.05, power=0.9, r=1:2, r0=1:2, p=0.75 , p0=0.5, u.shape="obf", l.shape="fixed", lfix=0, ufix=NULL, nstart=1, sample.size=TRUE, N=20){


  #require(mvtnorm) ## the function pmvnorm is required to evaluate multivariate normal probabilities

  ############################################################################################
  ## 'mesh' creates the points and respective weights to use in the outer quadrature integrals
  ## you provide one-dimensional set of points (x) and weights (w) (e.g. midpoint rule, gaussian quadrature, etc)
  ## and the number of dimensions (d) and it gives d-dimensional collection of points and weights
  ###############################################################################################

  mesh<-function(x,d,w=1/length(x)+x*0){
    n<-length(x)
    W<-X<-matrix(0,n^d,d)
    for (i in 1:d){
      X[,i]<-x
      W[,i]<-w
      x<-rep(x,rep(n,length(x)))
      w<-rep(w,rep(n,length(w)))
    }
    w<-exp(rowSums(log(W)))
    list(X=X,w=w)
  }

  prodsum<-function(x,l,u,r,r0,r0diff,J,K,Sigma){
    ########################################################################################################
    ## x is vector of dummy variables ( the t_j in gdt paper ), l and u are boundary vectors
    ########################################################################################################
    int<-prod(sapply(x,dnorm))
    L<-sqrt(r[1])/r0[1]*sqrt(r0diff[1])*x[1]+l[1]*sqrt(1+r[1]/r0[1])
    U<-sqrt(r[1])/r0[1]*sqrt(r0diff[1])*x[1]+u[1]*sqrt(1+r[1]/r0[1])
    insum<-pnorm(L)
    if(J>1){
      for (j in 2:J){
        L[j]<-sqrt(r[j])/r0[j]*(sqrt(r0diff[1:j])%*%x[1:j])+l[j]*sqrt(1+r[j]/r0[j])
        U[j]<-sqrt(r[j])/r0[j]*(sqrt(r0diff[1:j])%*%x[1:j])+u[j]*sqrt(1+r[j]/r0[j])
        insum<-insum+pmvnorm(lower=c(L[1:j-1],-Inf),upper=c(U[1:j-1],L[j]),sigma=Sigma[1:j,1:j])[1]
      }
    }
    int<-int*insum^K
    return(int)
  }

  #############################################################################################################
  ##  'typeI' performs the outer quadrature integral in type I error equation using 'mesh' and 'prodsum'
  ##  and calculates the difference with the nominal alpha.
  ##  The accuracy of the quadrature integral will depend on the choice of points and weights.
  ##  Here, the number of points in each dimension is an input, N.
  ##  The midpoint rule is used with a range of -6 to 6 in each dimension. 
  #############################################################################################################

  typeI<-function(C,alpha,N,r,r0,r0diff,J,K,Sigma,u.shape,l.shape,lfix=NULL,ufix=NULL){

    ########################################################################
    ## the form of the boundary constraints are determined as functions of C. 
    ######################################################################## 
    if(!is.function(u.shape)){
      if (u.shape=='obf'){
        u<-C*sqrt(r[J]/r)
      }
      else if (u.shape=='pocock'){
        u<-rep(C,J)
      }
      else if (u.shape=='fixed'){
        u<-c(rep(ufix,J-1),C)
      } 
      else if (u.shape=='triangular') {
        u<-C*(1+r/r[J])/sqrt(r) 
      }
    }else{
      u <- C*u.shape(J)
    }

    if(!is.function(l.shape)){
      if (l.shape=='obf'){
        l<- c(-C*sqrt(r[J]/r[1:(J-1)]),u[J])
      }
      else if (l.shape=='pocock'){
        l<-c(rep(-C,J-1),u[J])
      }
      else if (l.shape=='fixed'){
        l<-c(rep(lfix,J-1),u[J])
      } else if (l.shape=='triangular') {
        if(u.shape=="triangular"){
          l<--C*(1-3*r/r[J])/sqrt(r)
        }else{
          l<--C*(1-3*r/r[J])/sqrt(r)/(-1*(1-3)/sqrt(J))
        }
      }
    }else{
      l <- c(C*l.shape(J)[1:(J-1)],u[J])
    }


    mmp<-mesh((1:N-.5)/N*12-6,J,rep(12/N,N))
    evs<-apply(mmp$X,1,prodsum,l=l,u=u,r=r,r0=r0,r0diff=r0diff,J=J,K=K,Sigma=Sigma)
    truealpha<-1-mmp$w%*%evs
    return(truealpha-alpha)
  }

  ####################################################################################
  ## 'prodsum2' evaluates the integrand of Pi_1 according to  gdt paper:
  #####################################################################################

  prodsum2<-function(x,r,r0,l,u,K,delta,delta0,n,sig){
 
    int<-dnorm(x)  
    int<-int*pnorm(x+(delta-delta0)*sqrt(r[1]*n)/sig)^(K-1)*pnorm(sqrt(r0[1]/r[1])*(x+delta*sqrt(r[1]*n)/sig-u[1]*sqrt(1+r[1]/r0[1])))
    return(int)
  }

  ############################################################################################
  ## 'prodsum3' evaluates the integrand of Pi_j for j>1.
  ############################################################################################ 

  prodsum3 <- function(x,l,u,r,r0,r0diff,J,K,delta,delta0,n,sig,Sigma,SigmaJ){

    int<-prod(sapply(x,dnorm))
    L<-sqrt(r[1])/r0[1]*sqrt(r0diff[1])*x[1]+l[1]*sqrt(1+r[1]/r0[1])-delta0*sqrt(r[1]*n)/sig
    U<-sqrt(r[1])/r0[1]*sqrt(r0diff[1])*x[1]+u[1]*sqrt(1+r[1]/r0[1])-delta0*sqrt(r[1]*n)/sig
    insum<-pnorm(L)
    if (J>2){
      for (j in 2:(J-1)){
        L[j]<-sqrt(r[j])/r0[j]*(sqrt(r0diff[1:j])%*%x[1:j])+l[j]*sqrt(1+r[j]/r0[j])-delta0*sqrt(r[j]*n)/sig
        U[j]<-sqrt(r[j])/r0[j]*(sqrt(r0diff[1:j])%*%x[1:j])+u[j]*sqrt(1+r[j]/r0[j])-delta0*sqrt(r[j]*n)/sig
        insum<-insum+pmvnorm(lower=c(L[1:j-1],-Inf),upper=c(U[1:j-1],L[j]),sigma=Sigma[1:j,1:j])[1]
      }
    }
    U[J]<-x[J]+(delta-delta0)*sqrt(r[J]*n)/sig
    insum<-insum+pmvnorm(lower=c(L,-Inf),upper=U,sigma=Sigma)[1]
    int<-int*insum^(K-1)

    LJ<-sqrt(r[J]/(r[J]-r[1]))*(sqrt(r[1])/r0[1]*sqrt(r0diff[1])*x[1]+l[1]*sqrt(1+r[1]/r0[1])-delta*sqrt(r[1]*n)/sig-sqrt(r[1]/r[J])*x[J])
    UJ<-sqrt(r[J]/(r[J]-r[1]))*(sqrt(r[1])/r0[1]*sqrt(r0diff[1])*x[1]+u[1]*sqrt(1+r[1]/r0[1])-delta*sqrt(r[1]*n)/sig-sqrt(r[1]/r[J])*x[J])
    if (J>2){
      for (j in 2:(J-1)){
        LJ[j]<-sqrt(r[J]/(r[J]-r[j]))*(sqrt(r[j])/r0[j]*sqrt(r0diff[1:j])%*%x[1:j]+l[j]*sqrt(1+r[j]/r0[j])-delta*sqrt(r[j]*n)/sig-sqrt(r[j]/r[J])*x[J])
        UJ[j]<-sqrt(r[J]/(r[J]-r[j]))*(sqrt(r[j])/r0[j]*sqrt(r0diff[1:j])%*%x[1:j]+u[j]*sqrt(1+r[j]/r0[j])-delta*sqrt(r[j]*n)/sig-sqrt(r[j]/r[J])*x[J])
      }
    }
    int<-int*pmvnorm(lower=LJ,upper=UJ,sigma=SigmaJ)[1]
    int<-int*pnorm((r0[J]/sqrt(r[J])*(x[J]+delta*sqrt(r[J]*n)/sig-u[J]*sqrt(1+r[J]/r0[J]))-(sqrt(r0diff[1:(J-1)])%*%x[1:(J-1)]))/sqrt(r0diff[J]))
    return(int)

  }

  ######################################################################################################
  ##  'typeII' performs the outer quadrature integrals of Pi_j j=1,...,J  using 'mesh',
  ##  'prodsum2' and 'prodsum3' as well as summing the Pi_1,...,Pi_J and calculates the difference 
  ##  with the nominal power.
  ##  The accuracy of the quadrature integral again depends on the choice of points and weights.
  ##  Here, the number of points in each dimension is an input, N.
  ##  The midpoint rule is used with a range of -6 to 6 in each dimension. 
  ########################################################################################################

  typeII<-function(n,beta,l,u,N,r,r0,r0diff,J,K,delta,delta0,sig,Sigma){
  
    mmp<-mesh((1:N-.5)/N*12-6,1,rep(12/N,N))
    evs<-apply(mmp$X,1,prodsum2,r=r,r0=r0,l=l,u=u,K=K,delta=delta,delta0=delta0,n=n,sig=sig)
    pi<-mmp$w%*%evs

    if(J>1){
      for (j in 2:J){  
        if (j==2){
          A<-sqrt(r[j]/(r[j]-r[1:(j-1)]))
        }else{
          A<-diag(sqrt(r[j]/(r[j]-r[1:(j-1)])))
        }
        SigmaJ<-A%*%(Sigma[1:(j-1),1:(j-1)]-Sigma[1:(j-1),j]%*%t(Sigma[1:(j-1),j]))%*%A
        mmp<-mesh((1:N-.5)/N*12-6,j,rep(12/N,N))
        evs<-apply(mmp$X,1,prodsum3,l=l,u=u,r=r,r0=r0,r0diff=r0diff,J=j,K=K,delta=delta,delta0=delta0,n=n,sig=sig,Sigma=Sigma[1:j,1:j],SigmaJ=SigmaJ)
        pi<-pi+mmp$w%*%evs
      }
    }
    return(1-beta-pi)
  }

  ## checking input parameters
  if(K%%1>0 | J%%1>0){stop("K and J need to be integers.")}
  if(K < 1 | J < 1){stop("The number of stages and treatments must be at least 1.")}
  if(N<=3){stop("Number of points for integration by quadrature to small or negative.")}
  if(N>3 & N<=10){warning("Number of points for integration by quadrature is small which may result in inaccurate solutions.")}
  if(p<0 | p>1 | p0<0 | p0>1){stop("Treatment effect parameter not within 0 and 1.")}
  if(alpha<0 | alpha>1 | power<0 | power>1){stop("Error rate or power not between 0 and 1.")}
  if(p<p0){stop("Interesting treatment effect smaller than uninteresting effect.")}
  if(p0<0.5 ){warning("Uninteresting treatment effect less than 0.5 which implies that reductions in effect over placebo are interesting.")}
  if(length(r)!=length(r0)){stop('Different length of allocation ratios on control and experimental treatments')}
  if(length(r)!=J){stop('Length of allocation ratios does not match number of stages')}

  if(!is.function(u.shape)){
    if(!u.shape%in%c("pocock","obf","triangular","fixed")){stop("Upper boundary does not match the available options")}
    if(u.shape=="fixed" & is.null(ufix)){stop("ufix required when using a fixed upper boundary shape.")}
  }else{
    b <- u.shape(J)
    if(!all(sort(b,decreasing=TRUE)==b)){stop("Upper boundary shape is increasing")}
  }
  if(!is.function(l.shape)){
   if(!l.shape%in%c("pocock","obf","triangular","fixed")){stop("Lower boundary does not match the available options")}
   if(l.shape=="fixed" & is.null(lfix)){stop("lfix required when using a fixed lower boundary shape.")}
  }else{
    b <- l.shape(J)
    if(!all(sort(b,decreasing=FALSE)==b)){stop("Lower boundary shape is decreasing")}
  }
 
  ############################################################################
  ## Convert treatment effects into absolute effects with standard deviation 1:
  ############################################################################

  delta<-sqrt(2)*qnorm(p)
  delta0<-sqrt(2)*qnorm(p0)
  sig<-1
  
  ####################################################################
  ## Create the variance covariance matrix from allocation proportions:
  ####################################################################

  bottom<-matrix(r,J,J)
  top<-matrix(rep(r,rep(J,J)),J,J)
  top[upper.tri(top)]<-t(top)[upper.tri(top)]
  bottom[upper.tri(bottom)]<-t(bottom)[upper.tri(bottom)] 
  Sigma<-sqrt(top/bottom)

  ###############################################################################
  ## Create r0diff: the proportion of patients allocated to each particular stage
  ###############################################################################

  r0lag1<-c(0,r0[1:J-1])
  r0diff<-r0-r0lag1 


  ################################
  ## Find boundaries using 'typeI'
  ################################
  uJ<-NULL
  ## making sure that lfix is not larger then uJ
  try(uJ<-uniroot(typeI,c(qnorm(1-alpha)/2,5),alpha=alpha,N=N,r=r,r0=r0,r0diff=r0diff,J=J,K=K,Sigma=Sigma,u.shape=u.shape,l.shape=l.shape,lfix=lfix,ufix=ufix,tol=0.001)$root, silent=TRUE)
  if(is.null(uJ)){stop("Lower boundary (lfix) is too large.")}

  if(!is.function(u.shape)){
    if (u.shape=='obf'){
      u<-uJ*sqrt(r[J]/r)
    }
    else if (u.shape=='pocock'){
      u<-rep(uJ,J)
    }
    else if (u.shape=='fixed'){
      u<-c(rep(ufix,J-1),uJ)
    } else if (u.shape=='triangular') {
      u<-uJ*(1+r/r[J])/sqrt(r) 
    }
  }else{
    u <- uJ*u.shape(J)
  }

  if(!is.function(l.shape)){
    if (l.shape=='obf'){
      l<- c(-uJ*sqrt(r[J]/r[1:(J-1)]),u[J])
    }
    else if (l.shape=='pocock'){
      l<-c(rep(-uJ,J-1),u[J])
    }
    else if (l.shape=='fixed'){
      l<-c(rep(lfix,J-1),u[J])
    } else if (l.shape=='triangular') {
       if(u.shape=="triangular"){
          l<--uJ*(1-3*r/r[J])/sqrt(r)
       }else{
          l<--uJ*(1-3*r/r[J])/sqrt(r)/(-1*(1-3)/sqrt(J))
       }
    }
  }else{
    l <- c(uJ*l.shape(J)[1:(J-1)],u[J])
  }


  #########################################################
  ## Find alpha_star
  #########################################################
  alpha_star <- numeric(J)
  alpha_star[1] <- typeI(u[1], alpha = 0, N = N, r = r[1], r0 = r0[1], r0diff = r0diff[1], J = 1, K = K, Sigma = Sigma, u.shape = "fixed", l.shape = "fixed", lfix = NULL, ufix = NULL)
  if (J > 1){
      for (j in 2:J){
          alpha_star[j] <- typeI(u[j], alpha = 0, N = N, r = r[1:j], r0 = r0[1:j], r0diff = r0diff[1:j], J = j, K = K, Sigma = Sigma, u.shape = "fixed", l.shape = "fixed", lfix = l[1:(j - 1)], ufix = u[1:(j - 1)])
      }
  }


  
  #############################################################
  ##  Now find samplesize for arm 1 stage 1 (n)  using 'typeII'.
  ##  Sample sizes for all stages are then determined by
  ##  r*n and r0*n.
  #############################################################

  if(J==1 & p0==0.5) {
    if(r0>r){
      r <- r/r0
      r0 <- r0/r0
    }
    rho <- r / (r + r0)
    corr <- matrix(rho, K, K) + diag(1 - rho, K)
    quan <- qmvnorm(1-alpha, mean=rep(0, K), corr=corr)$quantile
    n <- ((quan + qnorm(power)) / (qnorm(p)*sqrt(2)))^2 *(1+1/r)

  }else{

    n<-nstart 
    ###################################################################################################
    ## program could be very slow starting at n=0, may want to start at more sensible lower bound on n
    ## unlike type I error, power equation does not neccessarily have unique solution n therefore search
    ## for smallest solution:
    ####################################################################################################
  
    pow<-0
    if(sample.size){
      while (pow==0){
        n<-n+1
        pow<-(typeII(n,beta=1-power,l=l,u=u,N=N,r=r,r0=r0,r0diff=r0diff,J=J,K=K,delta=delta,delta0=delta0,sig=sig,Sigma=Sigma)<0)
      }
    }else{
      n <- NULL
    }  
  }

  res <- NULL
  res$l <- l  
  res$u <- u
  res$n <- n

  res$rMat <- rbind(r0,matrix(r,ncol=J,nrow=K,byrow=TRUE)) ## allocation ratios
  res$N <- sum(ceiling(res$rMat[,J]*res$n)) ## maximum total sample sizeres$N <- K*r[J]*n+r0[J]*n ## maximum total sample size

  res$K <- K
  res$J <- J
  res$alpha <- alpha
  res$alpha_star <- alpha_star
  if(sample.size){
    res$power <- power
  }else{
    res$power <- NA
  }

  class(res)<-"MAMS"

  return(res)

}

####
## future functions
## - update boundaries (write separate function that takes a MAMS object)
## - correct for estimating variances using t-quantiles
## - expected sample size function
## - Implement code in C?
