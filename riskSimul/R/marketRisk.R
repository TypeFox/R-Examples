
########
#USED LIBRARY
########

#library("Runuran")
R2X<-function(tab,...){
  write.table(tab,"clipboard",sep="\t",row.names=F)
}

########
#NECESSARY FUNCTIONS
########

# OrthMat creates an orthogonal square matrix (basis) for given
# orthogonal directions (columns)

OrthMat<-function(Vbar  # D times Delta matrix that holds orthogonal directions
){
  D<-dim(Vbar)[1];Delta<-dim(Vbar)[2];
  V<-matrix(0,D,D);V[1:D,1:Delta]<-Vbar
  V[D,(Delta+1):D]<-1
  for(d in Delta:(D-1)){
    B<-V[(D-d):(D-1),1:d]
    b<-V[D,1:d]
    vbar<-solve(t(B),-b)
    V[(D-d):(D-1),d+1]<-vbar                 # orthogonal (d+1)st column determined
    V[,d+1]<-V[,d+1]/sqrt(sum(V[,d+1]^2))    # column is normalized
  }
  V                                          # return the final matrix
}

# OptAllocHeur is the optimal allocation heuristic that finds suboptimal
# solution for the minimization of J convex functions which are in the form of
# sum(a_i/pi_i) with a_i non-negative coefficients, the matrix A depends on the
# type of the error of simulation.

OptAllocHeur<-function(A,	# A is an I times J matrix
eps=1e-6				# critical error for convergence
){
  I<-dim(A)[1];J<-dim(A)[2]
  pi_j<-t(t(sqrt(A))/colSums(sqrt(A)))    # extreme points of the convex hull
  pi_c<-rowSums(pi_j)/J                   # current point: weight center
  delta<-1e16;eta<-1;pi_h<-pi_c           
  sol<-colSums(A/(pi_c+1e-16))            # current point all objective values
  j_c<-order(sol,decreasing=T)[1]         # improvement index
  omega_h<-omega_c<-sol[j_c]              # best and current objective value

  while((abs(delta)/omega_h) > eps){      
    pi_c<-(eta*pi_c+pi_j[,j_c])/(eta+1)   # move to candidate point
    sol<-colSums(A/(pi_c+1e-16))          # candidate point all objective values
    j_p<-order(sol,decreasing=T)[1]       # improvement index
    omega_p<-sol[j_p]                     # candidate objective value
    if(omega_p <= omega_h){			# if you have a better solution
      pi_h<-pi_c;omega_h<-omega_c         # update best solution
    }
    if(j_p != j_c){                       # if you have a new improvement index
      eta<-eta+1                          
      delta=omega_c-omega_p               
      omega_c<-omega_p;j_c<-j_p           # update current solution
    }
  }
  list(pi_h,omega_h,sol)                  # best solution, objective value, all solutions
}

# bvec calculates the stratified estimators for all thresholds.
# xbar is an array (J times I1 times I2) that holds conditional stratum means
# for J thresholds. For minimizing squared relative error, we calculate
# the simulation estimates as "b"

bvec<-function(xbar,    # J times I1 times I2 array holds cond. est. for each response 
p				# stratum probabilities (I1 times I2 matrix)
){
  d<-dim(xbar)[1]
  b<-sum(xbar[1,,]*p)^-2
  if(d>=2){
    for(j in 2:d){
      b[j]<-sum(xbar[j,,]*p)^-2
    }
  }
  b                     # a vector of length J, holds estimates for responses
}

# amat creates a J times I matrix (I=I1*I2) that holds conditional variances
# multiplied with squared stratum probabilities for each response. 

amat<-function(s2,	# J times I1 times I2 array holds cond. var. for each response
p				# stratum probabilities (I1 times I2 matrix)
){
  d<-dim(s2)[1]
  a<-as.vector(s2[1,,])*as.vector(p^2)
  if(d>=2){
    for(j in 2:d){
      a<-rbind(a,as.vector(s2[j,,])*as.vector(p^2))
    }
  }
  a				# a J times I matrix
}

#####xx###
#MODEL FUNCTIONS
########

# All simulation fuýnctions require a model and parameters object (list)
# this object has 6 variables

# 1) portfobj<-list("t" or "GH") indicates the type of marginal
# 2) portfobj[[2]]<-pmg or gen
#     either the parameters of marginal dist (as a matrix) 
#     or the generator object for uq function of Runuran
# 3) portfobj[[3]]<-nu dof of copula
# 4) portfobj[[4]]<-L lower triangular cholesky factorization of corr. matrix (D times D)
# 5) portfobj[[5]]<-c vector of scale factors of the portfolio (length D)
# 6) portfobj[[6]]<-w vector of investment fractions of protfolio (length D)

# ReturnCopula calculates the portfolio return for N replications of
# multinormal Z and gamma input Y. (Z: D times N, Y: length of N)

ReturnCopula<-function(Z, 	# multinormal input
Y,					# gamma (chi-squared) input
portfobj				# model parameters object
){
  D<-dim(Z)[1];N<-dim(Z)[2];                                    
  nu<-portfobj[[3]];L<-portfobj[[4]];c<-portfobj[[5]];w<-portfobj[[6]] # analysis of parameters
  T<-L%*%Z/matrix(sqrt(Y/nu),D,N,byrow=TRUE)                   # generating multi-T
  if(portfobj[[1]]=="t"){
    pmg<-portfobj[[2]]                                           
    marg<-qt(pt(T,nu),pmg[,3])*pmg[,2]+pmg[,1]                 # generate t marginals of t-copula
  }else if(portfobj[[1]]=="GH"){
    gen<-portfobj[[2]]
    marg<-matrix(0,D,N)
    for(i in 1:D){
      marg[i,]<-uq(gen[[i]],pt(T[i,],nu))                      # generate GH marginals of t-copula
    }
  }
  as.vector(t(w)%*%exp(c*marg))                                # weighted sum of log-scaled-marginals (portfolio return)
}

# TailLossProb simply returns 0 (or 1) if the return is less than (greater) than the treshold
# if there are J thresholds, returns a J times N matrix, each row for each threshold. 

TailLossProb<-function(return,	# vector of returns (size N)
threshold					# vector of thresholds (size J)
){
  if(length(threshold) > 1){
    return<-matrix(return,length(threshold),length(return),byrow=T) # replicate return for different thresholds
  }
  1*(return < threshold)                                            # a logical vector for loss level
}

# Excess simply returns 0 (or the loss) if the return is less than (greater) than the treshold
# if there are J thresholds, returns a J times N matrix, each row for each threshold. 

Excess<-function(return,	# vector of returns (size N)
threshold				# vector of thresholds (size J)
){
  if(length(threshold) > 1){
    return<-matrix(return,length(threshold),length(return),byrow=T) # replicate return for different thresholds
  }
  (1-return)*(return < threshold)                                   # loss levels for returns below threshold 
}

# evaluates return over a given direction, used for finding the root
# r for which we touch the rare event region for the first time

Touch<-function(r,	# the step size over given direction
dir,				# given direction
portfobj,			
threshold){   		
  Z<-t(t(r*dir));Y<-portfobj[[3]]                      
  ReturnCopula(Z,Y,portfobj)-threshold+1e-5
}

# Alg2 from Sak et al. 2010
# for given direction, finds optimal IS shift and scale parameters
# and the objective function 

Alg2<-function(dir,portfobj,threshold){
  dir<-dir/sqrt(sum(dir^2))                                  # normalize the direction 
  r0<-uniroot(Touch,c(-10^3,1e-5),dir,portfobj,threshold)$root # find the distance to the important region
  nu<-portfobj[[3]]                                            
  y0<-(nu-2)/(1+r0^2/nu)                                     
  z0<-r0*sqrt(y0/nu)*dir                                     
  objective<-(nu/2-1)*(log(y0)-1)                            # maximum IS density over the given direction
  list(c(z0,y0),objective) 
}

# Alg3 from Sak et al. 2010
# finds the optimal direction that yields minimum objective value
# returns optimal IS parameters

Alg3<-function(portfobj,threshold){
  L<-portfobj[[4]];c<-portfobj[[5]];w<-portfobj[[6]]    # seize the model parameters                                     # analysis of parameters
  dir<-as.vector(L%*%(c*w))				  # initial direction                                                            # intial direction
  dir<-dir/sqrt(sum(dir^2))                                                            # normalize the direction
  d<-length(dir)-1                                                                     # reduce dimension by one
  v<-dir[2:length(dir)]/dir[1]                    # reduce dimension by one                                             # scale over the first element
  f<-function(z){
    -Alg2(c(1,z),portfobj,threshold)[[2]]                                               
  }
  optdir<-optim(v,f,control=list(maxit=10000),lower=rep(0,d-1),method="L-BFGS-B")[[1]] # find optimal direction
  optdir<-c(1,optdir)                                                                  # add the first element                                                 
  optdir<-optdir/sqrt(sum(optdir^2))                                                   # normalize the optimal direction
  Alg2(optdir,portfobj,threshold)[[1]]                                                   # return z0 and y0 values
}


#################xx################
# NVCopula runs naive simulation for single threshold value
# Estimates both tail loss prob. and conditional excess

NVCopula<-function(N,	# size of the simulation
portfobj, 			# model parameter object
threshold,			# threshold
clock=F){			# if TRUE, returns also the execution time
  if(clock) start<-proc.time()                                                                    # timer of the algorithm
  D<-length(portfobj[[5]]);nu<-portfobj[[3]]                                                    
  if(portfobj[[1]]=="GH"){
    gen<-list()
    pmg <- portfobj[[2]]
    gen<-list()
    for(d in 1:D){
    	gen[[d]]<-pinvd.new(distr=udghyp(lambda=pmg[d,1], alpha=pmg[d,2], beta=pmg[d,3], delta=pmg[d,4], mu=pmg[d,5]))	
    }
    portfobj[[2]]<-gen
  }
  Z<-matrix(rnorm(D*N),D,N)                                                                        # multi normal input
  Y<-rgamma(N,shape=nu/2,scale=2)                                                                  # chi-squared input
  return<-ReturnCopula(Z,Y,portfobj)                                                                 # copula return
  yvec<-TailLossProb(return,threshold)                                                             # a logical vector of loss level
  xvec<-Excess(return,threshold)                                                                   # loss levels for returns below threshold 
  
  TLPROB<-EXCESS<-CONDEX<-rep(0,6)
  names(TLPROB)<-c("estimate","halfwidth","%95CILB","%95CIUB","est.var","rel.err")
  
  TLPROB[1]<-mean(yvec)                                                                            # probability estimator
  TLPROB[5]<-var(yvec)/N                                                                           # probability estimator variance
  EXCESS[1]<-mean(xvec)                                                                            # excess estimator 
  EXCESS[5]<-var(xvec)/N                                                                           # excess estimator variance
  SXY<-cov(yvec,xvec)/N
  CONDEX[1]<-EXCESS[1]/TLPROB[1]-(EXCESS[1]*TLPROB[5]/TLPROB[1]+SXY)/(TLPROB[1]^2)                 # conditional excess estimator
  CONDEX[5]<-TLPROB[5]*EXCESS[1]^2/TLPROB[1]^4-2*EXCESS[1]*SXY/(TLPROB[1]^3)+EXCESS[5]/TLPROB[1]^2 # conditional excess estimator variance
 
#  res<-rbind(TLPROB,EXCESS,CONDEX)
  res<-rbind(TLPROB,CONDEX)
  res[,2]<-sqrt(res[,5])*qnorm(0.975)                                                              # error bound
  res[,3]<-res[,1]-res[,2]                                                                         # conf. interval lower bound 
  res[,4]<-res[,1]+res[,2]                                                                         # conf. interval upper bound 
  res[,6]<-100*res[,2]/res[,1]                                                                     # relative error
  if(clock){
    finish<-proc.time()
    list(res,finish[3]-start[3])                                                                   # return results with timing
  }else{
    list(res)                                                                                      # return results without timing
  }
}

# searchx is the bisection method which find the threshold level for certain probabilities
# like p=0.05, 0.001, etc. It uses naive simulation.

searchx<-function(N,p,portfobj){
  a<-0.85;b<-1;
  fa<-NVCopula(N,portfobj,a)[[1]][1,1]-p
  fb<-NVCopula(N,portfobj,b)[[1]][1,1]-p
  err<-1
  while(err/p>0.05){
    ab<-(a+b)/2
    fab<-NVCopula(N,portfobj,ab)[[1]][1,1]-p
    err<-abs(fab)
    if(fab<0) a<-ab;fa<-fab
    if(fab>=0) b<-ab;fb<-fab
  }
  ab
}



#################xx#######
# NVCopulaMT runs naive simulation for multiple threshold values
# Estimates both tail loss prob. and conditional excess

NVCopulaMT<-function(N,	# size of the simulation
portfobj, 			# model parameter object
threshold,			# threshold vector
clock=F){			# if TRUE, returns also the execution time
  if(clock) start<-proc.time()                                                                     # timer of the algorithm
  D<-length(portfobj[[5]]);nu<-portfobj[[3]];J<-length(threshold)                                                   
  if(portfobj[[1]]=="GH"){
    pmg <- portfobj[[2]]
    gen<-list()
    for(d in 1:D){
    	gen[[d]]<-pinvd.new(distr=udghyp(lambda=pmg[d,1], alpha=pmg[d,2], beta=pmg[d,3], delta=pmg[d,4], mu=pmg[d,5]))	
    }
    portfobj[[2]]<-gen
  }
  Z<-matrix(rnorm(D*N),D,N)                                                                        # multi normal input
  Y<-rgamma(N,shape=nu/2,scale=2)                                                                  # chi-squared input
  return<-ReturnCopula(Z,Y,portfobj)                                                                 # copula return
  ymat<-TailLossProb(return,threshold)                                                             # a logical matrix of loss level
  xmat<-Excess(return,threshold)                                                                   # loss levels for returns below threshold 
  
  TLPROB<-EXCESS<-CONDEX<-matrix(0,J,6)
  colnames(TLPROB)<-c("estimate","halfwidth","%95CILB","%95CIUB","est.var","rel.err")
  for(j in 1:J){
    TLPROB[j,1]<-mean(ymat[j,])
    TLPROB[j,5]<-var(ymat[j,])/N                                                                           # probability estimator variance
    EXCESS[j,1]<-mean(xmat[j,])                                                                            # excess estimator 
    EXCESS[j,5]<-var(xmat[j,])/N                                                                           # excess estimator variance
    SXY<-cov(ymat[j,],xmat[j,])/N
    CONDEX[j,1]<-EXCESS[j,1]/TLPROB[j,1]-(EXCESS[j,1]*TLPROB[j,5]/TLPROB[j,1]+SXY)/(TLPROB[j,1]^2)                 # conditional excess estimator
    CONDEX[j,5]<-TLPROB[j,5]*EXCESS[j,1]^2/TLPROB[j,1]^4-2*EXCESS[j,1]*SXY/(TLPROB[j,1]^3)+EXCESS[j,5]/TLPROB[j,1]^2 # conditional excess estimator variance
  }

  res<-rbind(TLPROB,EXCESS,CONDEX)
  res[,2]<-sqrt(res[,5])*qnorm(0.975)                                                              # error bound
  res[,3]<-res[,1]-res[,2]                                                                         # conf. interval lower bound 
  res[,4]<-res[,1]+res[,2]                                                                         # conf. interval upper bound 
  res[,6]<-100*res[,2]/res[,1]                                                                     # relative error
  if(clock){
    finish<-proc.time()
    result<-list(res[1:J,],res[J+(1:J),],res[2*J+(1:J),],threshold,finish[3]-start[3])                                                                   # return results with timing
    names(result)<-c("TLPROB","EXCESS","CONDEX","THRESHOLDS","TIME")
  }else{
#    result<-list(res[1:J,],res[J+(1:J),],res[2*J+(1:J),],threshold)                                                                                      # return results without timing
#    names(result)<-c("TLPROB","EXCESS","CONDEX","THRESHOLDS")
    result<-list(res[1:J,],res[2*J+(1:J),],threshold)                                                                                      # return results without timing
    names(result)<-c("TLPROB","CONDEX","THRESHOLDS")
  }
  result
}


#############xx#######
# SISCopula returns stratified estimators of tail loss prob. and conditional excess for a single threshold value

SISCopula<-function(N,	# a vector of sample allocations through iterations
stratasize,			# vector of length two, holds number of strata for each input e.g. c(22,22)
portfobj,          		# model parameters object
threshold,			# theshold value
CEopt=F,			# if TRUE, minimizes the variance of conditional excess estimate
clock=F
){
  if(clock) start<-proc.time()                                                                                     # timer of the algorithm
#### DIRECTION OPTIMIZATION
  D<-length(portfobj[[5]]);nu<-portfobj[[3]]
  if(portfobj[[1]]=="GH"){															 # creates generator objects for GH marginals
    pmg <- portfobj[[2]]
    gen<-list()
    for(d in 1:D){
    	gen[[d]]<-pinvd.new(distr=udghyp(lambda=pmg[d,1], alpha=pmg[d,2], beta=pmg[d,3], delta=pmg[d,4], mu=pmg[d,5]))	
    }
    portfobj[[2]]<-gen
  }
  muy0<-Alg3(portfobj,threshold)                                                                       			 # finds optimal IS parameters and obj val.	
  theta<-muy0[D+1]/(nu/2-1);mu<-muy0[1:D]                                                                          # optimal IS parameters
  v<-t(t(mu))                           
  mu<-sqrt(sum(mu^2))                                                                                              # IS shift (after linear transformation)
  v<-v/mu                                                                                                          # stratification direction                                                          
  V<-OrthMat(v)                                                                                                    # linear transformation matrix                                                                    
#### STRATA CONSTRUCTION            
  I1<-stratasize[1];I2<-stratasize[2]                                                                              # strata sizes
  strata1<-0:I1/I1;strata2<-0:I2/I2                                                                                # strata boundaries
  I<-I1*I2                                                                                                         # total strata size
#### INITIALIZATION
  p<-1/I;																		 # stratum probabilities
  xbar<-ybar<-s2x<-s2y<-sxy<-matrix(0,I1,I2)												 # conditional mean and var. arrays
# x is for the numerator, y is for the denominator of cond. excess estimate (y is for tail loss prob.)
#### FIRST ITERATION
  m<-p*N[1]																		 # proportional allocation in the first iter.
  fmcumsum<-floor(cumsum(as.vector(m)))
  M<-matrix(c(fmcumsum[1],diff(fmcumsum)),I1,I2)											 # force allocations to be integer (M is I1*I2 matrix)
  M<-(M<10)*10+(M>=10)*M															 # force at least 10 in each stratum

  StU1<-rep(strata1[-(I1+1)],rowSums(M))+runif(sum(M))*rep(diff(strata1),rowSums(M))                               # strf. unif. for multi normal input
  Z<-matrix(c(qnorm(StU1,mu),rnorm(sum(M)*(D-1))),D,sum(M),byrow=TRUE)                                             # inversion to normal 
  StU2<-rep(rep(strata2[-(I2+1)],I1),as.vector(t(M)))+runif(sum(M))*rep(rep(diff(strata2),I1),as.vector(t(M)))     # double strf. unif. for gamma input
  Y<-qgamma(StU2,shape=nu/2,scale=theta)                                                                           # inversion to gamma 
  wIS<-exp((0.5*mu-Z[1,])*mu+(2-theta)*Y/(2*theta))*(theta/2)^(nu/2)                                               # IS weight (ratio) 
  return<-ReturnCopula(V%*%Z,Y,portfobj)                                                                             # copula return
  yvec<-wIS*TailLossProb(return,threshold)                                                                         # a logical vector of loss level
  xvec<-wIS*Excess(return,threshold)                                                                               # loss levels for returns below threshold                                       

  Mcumsum<-cumsum(t(M))
  sampley<-list(list(yvec[1:Mcumsum[1]]))													 # sample set for denominator in stratum one
  samplex<-list(list(xvec[1:Mcumsum[1]]))													 # sample set for numerator in stratum one	
  ybar[1,1]<-mean(sampley[[1]][[1]]);s2y[1,1]<-var(sampley[[1]][[1]])                                              # xbar ybar holds conditional estimates
  xbar[1,1]<-mean(samplex[[1]][[1]]);s2x[1,1]<-var(samplex[[1]][[1]])                                              # s2x s2y holds conditional variances
  sxy[1,1]<-cov(sampley[[1]][[1]],samplex[[1]][[1]])                                                               # sxy holds conditional covariances
  for(j in 2:I2){																	 # do the same for other strata
    sampley[[1]][[j]]<-yvec[(Mcumsum[j-1]+1):Mcumsum[j]]
    samplex[[1]][[j]]<-xvec[(Mcumsum[j-1]+1):Mcumsum[j]]
    ybar[1,j]<-mean(sampley[[1]][[j]]);s2y[1,j]<-var(sampley[[1]][[j]])
    xbar[1,j]<-mean(samplex[[1]][[j]]);s2x[1,j]<-var(samplex[[1]][[j]])
    sxy[1,j]<-cov(sampley[[1]][[j]],samplex[[1]][[j]])
  }
  for(i in 2:I1){																	 # do the same for other strata
    sampley[[i]]<-samplex[[i]]<-list()
    for(j in 1:I2){
      sampley[[i]][[j]]<-yvec[(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]]
      samplex[[i]][[j]]<-xvec[(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]]
      ybar[i,j]<-mean(sampley[[i]][[j]]);s2y[i,j]<-var(sampley[[i]][[j]])
      xbar[i,j]<-mean(samplex[[i]][[j]]);s2x[i,j]<-var(samplex[[i]][[j]])
      sxy[i,j]<-cov(sampley[[i]][[j]],samplex[[i]][[j]])
    }
  }
  n<-M
#### ADAPTIVE ITERATIONS
  if(length(N)>=2){																 
    for(k in 2:length(N)){															 # in adaptive iteration k of AOA algorithm
      if(CEopt==F){																 # opt. alloc. for tail loss prob.
        sy<-sqrt(s2y)
        m<-N[k]*sy/sum(sy)
      }else{																	 # opt. alloc. for conditional excess
        x<-p*sum(xbar);y<-p*sum(ybar)
        sr<-sqrt(x^2*s2y/y^4 - 2*x*sxy/y^3 + s2x/y^2)
        m<-N[k]*sr/sum(sr)
      }
      fmcumsum<-floor(cumsum(as.vector(m)))
      M<-matrix(c(fmcumsum[1],diff(fmcumsum)),I1,I2)
      M<-(M<10)*10+(M>=10)*M

      StU1<-rep(strata1[-(I1+1)],rowSums(M))+runif(sum(M))*rep(diff(strata1),rowSums(M))                           # strf. unif. for multi normal input
      Z<-matrix(c(qnorm(StU1,mu),rnorm(sum(M)*(D-1))),D,sum(M),byrow=TRUE)                                         # inversion to normal ?? pinv.new
      StU2<-rep(rep(strata2[-(I2+1)],I1),as.vector(t(M)))+runif(sum(M))*rep(rep(diff(strata2),I1),as.vector(t(M))) # double strf. unif. for gamma input
      Y<-qgamma(StU2,shape=nu/2,scale=theta)                                                                       # inversion to gamma ?? pinv.new
      wIS<-exp((0.5*mu-Z[1,])*mu+(2-theta)*Y/(2*theta))*(theta/2)^(nu/2)                                           # IS weight (ratio) ??
      return<-ReturnCopula(V%*%Z,Y,portfobj)                                                                         # copula return
      yvec<-wIS*TailLossProb(return,threshold)                                                                     # a logical vector of loss level
      xvec<-wIS*Excess(return,threshold)                                                                           # loss levels for returns below threshold

      Mcumsum<-cumsum(t(M))
      sampley[[1]][[1]]<-c(sampley[[1]][[1]],yvec[1:Mcumsum[1]])
      samplex[[1]][[1]]<-c(samplex[[1]][[1]],xvec[1:Mcumsum[1]])
      ybar[1,1]<-mean(sampley[[1]][[1]]);s2y[1,1]<-var(sampley[[1]][[1]])
      xbar[1,1]<-mean(samplex[[1]][[1]]);s2x[1,1]<-var(samplex[[1]][[1]])
      sxy[1,1]<-cov(sampley[[1]][[1]],samplex[[1]][[1]])
      for(j in 2:I2){
        sampley[[1]][[j]]<-c(sampley[[1]][[j]],yvec[(Mcumsum[j-1]+1):Mcumsum[j]])
        samplex[[1]][[j]]<-c(samplex[[1]][[j]],xvec[(Mcumsum[j-1]+1):Mcumsum[j]])
        ybar[1,j]<-mean(sampley[[1]][[j]]);s2y[1,j]<-var(sampley[[1]][[j]])
        xbar[1,j]<-mean(samplex[[1]][[j]]);s2x[1,j]<-var(samplex[[1]][[j]])
        sxy[1,j]<-cov(sampley[[1]][[j]],samplex[[1]][[j]])
      }
      for(i in 2:I1){
        for(j in 1:I2){
          sampley[[i]][[j]]<-c(sampley[[i]][[j]],yvec[(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]])
          samplex[[i]][[j]]<-c(samplex[[i]][[j]],xvec[(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]])
          ybar[i,j]<-mean(sampley[[i]][[j]]);s2y[i,j]<-var(sampley[[i]][[j]])
          xbar[i,j]<-mean(samplex[[i]][[j]]);s2x[i,j]<-var(samplex[[i]][[j]])
          sxy[i,j]<-cov(sampley[[i]][[j]],samplex[[i]][[j]])
        }
      }
      n<-n+M
    }
  }
#### OUTPUT
  TLPROB<-EXCESS<-CONDEX<-rep(0,6)
  names(TLPROB)<-c("estimate","halfwidth","%95CILB","%95CIUB","est.var","rel.err")
  
  TLPROB[1]<-p*sum(ybar)                                                                                           # probability estimator
  TLPROB[5]<-p^2*sum(s2y/n)                                                                                        # probability estimator variance
  EXCESS[1]<-p*sum(xbar)                                                                                           # excess estimator 
  EXCESS[5]<-p^2*sum(s2x/n)                                                                                        # excess estimator variance
  SXY<-p^2*sum(sxy/n)
  CONDEX[1]<-EXCESS[1]/TLPROB[1]-(EXCESS[1]*TLPROB[5]/TLPROB[1]+SXY)/(TLPROB[1]^2)                                 # conditional excess estimator
  CONDEX[5]<-TLPROB[5]*EXCESS[1]^2/TLPROB[1]^4-2*EXCESS[1]*SXY/(TLPROB[1]^3)+EXCESS[5]/TLPROB[1]^2                 # conditional excess estimator variance
 
#  res<-rbind(TLPROB,EXCESS,CONDEX)
  res<-rbind(TLPROB,CONDEX)
  res[,2]<-sqrt(res[,5])*qnorm(0.975)                                                                              # error bound
  res[,3]<-res[,1]-res[,2]                                                                                         # conf. interval lower bound 
  res[,4]<-res[,1]+res[,2]                                                                                         # conf. interval upper bound 
  res[,6]<-100*res[,2]/res[,1]                                                                                     # relative error
  if(clock){
    finish<-proc.time()
    list(res,finish[3]-start[3])                                                                                   # return results with timing
  }else{
    res                                                                                                      # return results without timing
  }
}


###########xx####################################################################################################
# SISCopulaMT returns stratified estimators of tail loss prob. and conditional excess for multiple threshold values
# Choose either to minimize the overall error of tail loss prob. estimates
# or the overall error of conditional excess estimates
# The IS parameter is selected for an intermediate threshold value

SISCopulaMT<-function(N,			# a vector that holds allocations through iterations c(1,4,5)*10^5
stratasize,						# a vector of length two, holds strata sizes for each random input c(22,22)
portfobj,						# model paramters object
threshold,						# vector of thresholds (possibly ordered)
beta=0.75,						# coefficient of max. threshold value to choose an intermediate threshold
mintype,						# if 0: MSE, -1: MSR, -2: MAXE, -3: MAXR
							# if j (positive integer), then minimize the var. of the j-th estimate
CEopt=F,						# if TRUE, minimize the overall error of CE estimates
clock=F						# if true, return execution time
){
  if(clock) start<-proc.time()                                                                                     # timer of the algorithm
#### DIRECTION OPTIMIZATION
  D<-length(portfobj[[5]]);nu<-portfobj[[3]];J<-length(threshold)
  if(portfobj[[1]]=="GH"){
    pmg <- portfobj[[2]]
    gen<-list()
    for(d in 1:D){
    	gen[[d]]<-pinvd.new(distr=udghyp(lambda=pmg[d,1], alpha=pmg[d,2], beta=pmg[d,3], delta=pmg[d,4], mu=pmg[d,5]))	
    }
    portfobj[[2]]<-gen
  }
  thresholdstar<-max(threshold)*beta+min(threshold)*(1-beta)									 # optimal IS parameters for intermediate threshold
  muy0<-Alg3(portfobj,thresholdstar)                                                                       		 # finds optimal IS parameters and obj val.
  theta<-muy0[D+1]/(nu/2-1);mu<-muy0[1:D]                                                                          # optimal IS parameters
  v<-t(t(mu))                           
  mu<-sqrt(sum(mu^2))                                                                                              # IS shift (after transformation)
  v<-v/mu                                                                                                          # stratification direction                                                          
  V<-OrthMat(v)                                                                                                    # linear transformation matrix                                                                    
#### STRATA CONSTRUCTION            
  I1<-stratasize[1];I2<-stratasize[2]                                                                              # strata sizes
  strata1<-0:I1/I1;strata2<-0:I2/I2                                                                                # strata boundaries
  I<-I1*I2                                                                                                         # total strata size
#### INITIALIZATION
  p<-1/I;																		 # stratum probabilities
  xbar<-ybar<-s2x<-s2y<-sxy<-array(0,c(J,I1,I2))											 # arrays hold cond. mean, var and cov. for the
																			 # the numerator and the denominator for each response  
#### FIRST ITERATION
  m<-p*N[1]
  fmcumsum<-floor(cumsum(as.vector(m)))
  M<-matrix(c(fmcumsum[1],diff(fmcumsum)),I1,I2)
  M<-(M<10)*10+(M>=10)*M

  StU1<-rep(strata1[-(I1+1)],rowSums(M))+runif(sum(M))*rep(diff(strata1),rowSums(M))                               # strf. unif. for multi normal input
  Z<-matrix(c(qnorm(StU1,mu),rnorm(sum(M)*(D-1))),D,sum(M),byrow=TRUE)                                             # inversion to normal 
  StU2<-rep(rep(strata2[-(I2+1)],I1),as.vector(t(M)))+runif(sum(M))*rep(rep(diff(strata2),I1),as.vector(t(M)))     # double strf. unif. for gamma input
  Y<-qgamma(StU2,shape=nu/2,scale=theta)                                                                           # inversion to gamma 
  wIS<-exp((0.5*mu-Z[1,])*mu+(2-theta)*Y/(2*theta))*(theta/2)^(nu/2)                                               # IS weight (ratio) 
  return<-ReturnCopula(V%*%Z,Y,portfobj)                                                                             # copula return
  ymat<-t(wIS*t(TailLossProb(return,threshold)))                                                                   # a logical matrix (J*sample size) of loss level (denominator)
  xmat<-t(wIS*t(Excess(return,threshold)))												 # a matrix that holds excess values (numerator)

  Mcumsum<-cumsum(t(M))
  samplex<-sampley<-list()
  for(k in 1:J){																	 # for each threshold
    sampley[[k]]<-list(list(ymat[k,1:Mcumsum[1]]))											 # sample set of stratum 1 for the denominator
    samplex[[k]]<-list(list(xmat[k,1:Mcumsum[1]]))											 # sample set of stratum 1 for the numerator
    ybar[k,1,1]<-mean(sampley[[k]][[1]][[1]]);s2y[k,1,1]<-var(sampley[[k]][[1]][[1]])                              # xbar ybar holds conditional means
    xbar[k,1,1]<-mean(samplex[[k]][[1]][[1]]);s2x[k,1,1]<-var(samplex[[k]][[1]][[1]])					 # s2x and s2y holds conditional variances
    sxy[k,1,1]<-cov(sampley[[k]][[1]][[1]],samplex[[k]][[1]][[1]])								 # sxy holds cond. covariances
    for(j in 2:I2){																 # repeat for other strata	
      sampley[[k]][[1]][[j]]<-ymat[k,(Mcumsum[j-1]+1):Mcumsum[j]]
      samplex[[k]][[1]][[j]]<-xmat[k,(Mcumsum[j-1]+1):Mcumsum[j]]
      ybar[k,1,j]<-mean(sampley[[k]][[1]][[j]]);s2y[k,1,j]<-var(sampley[[k]][[1]][[j]])
      xbar[k,1,j]<-mean(samplex[[k]][[1]][[j]]);s2x[k,1,j]<-var(samplex[[k]][[1]][[j]])
      sxy[k,1,j]<-cov(sampley[[k]][[1]][[j]],samplex[[k]][[1]][[j]])
    }
    for(i in 2:I1){																 # repeat for other strata
      sampley[[k]][[i]]<-samplex[[k]][[i]]<-list()
      for(j in 1:I2){
        sampley[[k]][[i]][[j]]<-ymat[k,(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]]
        samplex[[k]][[i]][[j]]<-xmat[k,(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]]
        ybar[k,i,j]<-mean(sampley[[k]][[i]][[j]]);s2y[k,i,j]<-var(sampley[[k]][[i]][[j]])
        xbar[k,i,j]<-mean(samplex[[k]][[i]][[j]]);s2x[k,i,j]<-var(samplex[[k]][[i]][[j]])
        sxy[k,i,j]<-cov(sampley[[k]][[i]][[j]],samplex[[k]][[i]][[j]])
      }
    }
  }
  n<-M
#### ADAPTIVE ITERATIONS
  if(length(N)>=2){																 # adaptive iterations of the AOA algorithm										
    for(l in 2:length(N)){	
      if(CEopt==F){															# if min. for tail loss prob.
        b<-bvec(ybar,p)															 # vector of tlp estimates (size J)
        a<-amat(s2y,p)															 # matrix of cond.var. of tlp (J times I)	
      }else{																# if min. for cond. excess
        b<-bvec(xbar,p)/bvec(ybar,p)												 # vector of cond.excess estimates (size J)	
        s2r<-array(0,c(J,I1,I2))													 # creates a vector of conditional variances	
        for(k in 1:J){
          x<-p*sum(xbar[k,,]);y<-p*sum(ybar[k,,]);
          s2r[k,,]<-x^2*s2y[k,,]/y^4 - 2*x*sxy[k,,]/y^3 + s2x[k,,]/y^2							 # allocation fractions should be proportional to s2r
        }
        a<-amat(s2r,p)															 # matrix of s2r values for cond.exc. (J times I)
      }
      if(mintype == 0){															 # minimize MSE
        m<-N[l]*sqrt(colSums(a))/sum(sqrt(colSums(a)))
      }else if(mintype== -1){														 # minimize MSR
        m<-N[l]*sqrt(colSums(a*b))/sum(sqrt(colSums(a*b)))
      }else if(mintype== -2){														 # minimize MAXE
        m<-N[l]*matrix(OptAllocHeur(t(a),1e-6)[[1]],I1,I2)
      }else if(mintype== -3){														 # minimize MAXR
        m<-N[l]*matrix(OptAllocHeur(t(b*a),1e-6)[[1]],I1,I2)
      }else{																 # minimize the j-th estimate 		
        m<-N[l]*sqrt(s2r[mintype,,])/sum(sqrt(s2r[mintype,,]))
      }
      fmcumsum<-floor(cumsum(as.vector(m)))
      M<-matrix(c(fmcumsum[1],diff(fmcumsum)),I1,I2)
      M<-(M<10)*10+(M>=10)*M

      StU1<-rep(strata1[-(I1+1)],rowSums(M))+runif(sum(M))*rep(diff(strata1),rowSums(M))                           # strf. unif. for multi normal input
      Z<-matrix(c(qnorm(StU1,mu),rnorm(sum(M)*(D-1))),D,sum(M),byrow=TRUE)                                         # inversion to normal 
      StU2<-rep(rep(strata2[-(I2+1)],I1),as.vector(t(M)))+runif(sum(M))*rep(rep(diff(strata2),I1),as.vector(t(M))) # double strf. unif. for gamma input
      Y<-qgamma(StU2,shape=nu/2,scale=theta)                                                                       # inversion to gamma 
      wIS<-exp((0.5*mu-Z[1,])*mu+(2-theta)*Y/(2*theta))*(theta/2)^(nu/2)                                           # IS weight (ratio)
      return<-ReturnCopula(V%*%Z,Y,portfobj)                                                                         # copula return
      ymat<-t(wIS*t(TailLossProb(return,threshold)))                                                               # a logical vector of loss level
      xmat<-t(wIS*t(Excess(return,threshold)))                                                                     # loss levels for returns below threshold

      Mcumsum<-cumsum(t(M))
      for(k in 1:J){																 # combine sample sets with the current sample
        sampley[[k]][[1]][[1]]<-c(sampley[[k]][[1]][[1]],ymat[k,1:Mcumsum[1]])
        samplex[[k]][[1]][[1]]<-c(samplex[[k]][[1]][[1]],xmat[k,1:Mcumsum[1]])
        ybar[k,1,1]<-mean(sampley[[k]][[1]][[1]]);s2y[k,1,1]<-var(sampley[[k]][[1]][[1]])
        xbar[k,1,1]<-mean(samplex[[k]][[1]][[1]]);s2x[k,1,1]<-var(samplex[[k]][[1]][[1]])
        sxy[k,1,1]<-cov(sampley[[k]][[1]][[1]],samplex[[k]][[1]][[1]])
        for(j in 2:I2){
          sampley[[k]][[1]][[j]]<-c(sampley[[k]][[1]][[j]],ymat[k,(Mcumsum[j-1]+1):Mcumsum[j]])
          samplex[[k]][[1]][[j]]<-c(samplex[[k]][[1]][[j]],xmat[k,(Mcumsum[j-1]+1):Mcumsum[j]])
          ybar[k,1,j]<-mean(sampley[[k]][[1]][[j]]);s2y[k,1,j]<-var(sampley[[k]][[1]][[j]])
          xbar[k,1,j]<-mean(samplex[[k]][[1]][[j]]);s2x[k,1,j]<-var(samplex[[k]][[1]][[j]])
          sxy[k,1,j]<-cov(sampley[[k]][[1]][[j]],samplex[[k]][[1]][[j]])
        }
        for(i in 2:I1){
          for(j in 1:I2){
            sampley[[k]][[i]][[j]]<-c(sampley[[k]][[i]][[j]],ymat[k,(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]])
            samplex[[k]][[i]][[j]]<-c(samplex[[k]][[i]][[j]],xmat[k,(Mcumsum[(i-1)*I2+j-1]+1):Mcumsum[(i-1)*I2+j]])
            ybar[k,i,j]<-mean(sampley[[k]][[i]][[j]]);s2y[k,i,j]<-var(sampley[[k]][[i]][[j]])
            xbar[k,i,j]<-mean(samplex[[k]][[i]][[j]]);s2x[k,i,j]<-var(samplex[[k]][[i]][[j]])
            sxy[k,i,j]<-cov(sampley[[k]][[i]][[j]],samplex[[k]][[i]][[j]])
          }
        }
      }
      n<-n+M
    }
  }
#### OUTPUT
  TLPROB<-EXCESS<-CONDEX<-matrix(0,J,6)
  for(j in 1:J){
    TLPROB[j,1]<-p*sum(ybar[j,,])
    TLPROB[j,5]<-p^2*sum(s2y[j,,]/n)                                                                         # probability estimator variance
    EXCESS[j,1]<-p*sum(xbar[j,,])                                                                            # excess estimator 
    EXCESS[j,5]<-p^2*sum(s2x[j,,]/n)                                                                         # excess estimator variance
    SXY<-p^2*sum(sxy[j,,]/n)
    CONDEX[j,1]<-EXCESS[j,1]/TLPROB[j,1]-(EXCESS[j,1]*TLPROB[j,5]/TLPROB[j,1]+SXY)/(TLPROB[j,1]^2)           # conditional excess estimator
    CONDEX[j,5]<-TLPROB[j,5]*EXCESS[j,1]^2/TLPROB[j,1]^4-2*EXCESS[j,1]*SXY/(TLPROB[j,1]^3)+EXCESS[j,5]/TLPROB[j,1]^2 # conditional excess estimator variance
  }
  colnames(TLPROB)<-c("estimate","halfwidth","%95CILB","%95CIUB","est.var","rel.err")
  res<-rbind(TLPROB,EXCESS,CONDEX)
  res[,2]<-sqrt(res[,5])*qnorm(0.975)                                                              # error bound
  res[,3]<-res[,1]-res[,2]                                                                         # conf. interval lower bound 
  res[,4]<-res[,1]+res[,2]                                                                         # conf. interval upper bound 
  res[,6]<-100*res[,2]/res[,1]                                                                     # relative error
  if(clock){
    finish<-proc.time()
    result<-list(res[1:J,],res[J+(1:J),],res[2*J+(1:J),],threshold,finish[3]-start[3])             # return results with timing
    names(result)<-c("TLPROB","EXCESS","CONDEX","THRESHOLDS","TIME")
  }else{
#    result<-list(res[1:J,],res[J+(1:J),],res[2*J+(1:J),],threshold)                                # return results without timing
#    names(result)<-c("TLPROB","EXCESS","CONDEX","THRESHOLDS")
    result<-list(res[1:J,],res[2*J+(1:J),],threshold)                                # return results without timing
    names(result)<-c("TLPROB","CONDEX","THRESHOLDS")
  }
  result
}


#############################
#
# wrapper functions that will be exported

#################################
NVTCopula <- function(n=10^5,  # total sample size		
portfobj,  						# model paramters object
threshold =c(0.95,0.9) 		# vector of thresholds (possibly ordered)
){
 if(length(threshold)>1){
   res<- NVCopulaMT(N=n,portfobj=portfobj,threshold=threshold,clock=F)
 }else{
   res<- NVCopula(N=n,portfobj=portfobj,threshold=threshold,clock=F)
 }
return(res)
}


###########xx####################################################################################################
# SISCopulaMT returns stratified estimators of tail loss prob. and conditional excess for multiple threshold values
# Choose either to minimize the overall error of tail loss prob. estimates
# or the overall error of conditional excess estimates
# The IS parameter is selected for an intermediate threshold value


SISTCopula <- function(n=10^5,  # total sample size		
npilot=c(10^4,2*10^4),			# a vector that holds the size of the pilot iterations c(10^4,2*10^4)
portfobj,  						# object of portfolio parameters
threshold=c(0.95,0.9),			# vector of thresholds (possibly ordered)
stratasize=c(22,22),			# a vector of length two, holds strata sizes for each random input c(22,22)
CEopt=FALSE,					# if TRUE, minimize the overall error of CE estimates
    # next two arguments only for length(threshold) > 1
beta=0.75,						# coefficient of max. threshold value to choose an intermediate threshold
mintype=-1 						# if 0: MSE, -1: MSR, -2: MAXE, -3: MAXR
							    # if j (positive integer), then minimize the variance of the estimate for the j-th threshold 
){
 if(length(threshold)>1){
   res<- SISCopulaMT(N=c(npilot,n-sum(npilot)),stratasize=stratasize,portfobj=portfobj,threshold=threshold,beta=beta,mintype=mintype,CEopt=CEopt,clock=F)
 }else{
   res<- SISCopula(N=c(npilot,n-sum(npilot)),stratasize=stratasize,portfobj=portfobj,threshold=threshold,CEopt=CEopt,clock=F) 
 }
return(res)
}

new.portfobj <- function(nu,R,typemg="GH",parmg,c=rep(1,dim(R)[1]),w=c/sum(c)){
 return(list(typemg=typemg,parmg=parmg,nu=nu,L=t(chol(R)),c=c,w=w))
}
