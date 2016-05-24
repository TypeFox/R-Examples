# @n stands for number of agents
# @m stands for number of alternatives
# @distribution can be Normal or Exponential

# The Code Contains 4 different parts
# 1) Data Generation
# 2) Model Likelihood
# 3) Estimation Methods
# 4) Inference

################################################################################
# Data Generation
################################################################################

#' Parameter Generation for a RUM model
#' 
#' Exponential models mean parameters are drawn from a uniform distribution
#' Normal models, mean and standard devaition parameters are drawn from a standard unifrom
#' 
#' @param m number of sets of parameters to be drawn
#' @param distribution either 'normal' or 'exponential'
#' @return a list of RUM parameters
#' @export
#' @examples
#' Generate.RUM.Parameters(10, "normal")
#' Generate.RUM.Parameters(10, "exponential")

Generate.RUM.Parameters= function(m, distribution) {
  if(distribution=='normal'){parameter=list(m=m,Mean=runif(m,0,1),SD=runif(m,0,1))}
  else if(distribution=='exponential'){
    unscaled  <-  runif(m,0,1)
    parameter <- list(m = m, Mean=unscaled/sum(unscaled))
  }
  else stop(paste0("Distribution name \"", distribution, "\" not recognized"))
  parameter
}

#' Generate observation of ranks given parameters
#' 
#' Given a list of parameters (generated via the Generate RUM Parameters function),
#' generate random utilities from these models and then return their ranks
#' 
#' @param Params inference object from an Estimation function, or parameters object from a generate function
#' @param m number of alternatives
#' @param n number of agents
#' @param distribution can be either 'normal' or 'exponential'
#' @return a matrix of observed rankings
#' @export
#' @examples
#' Params = Generate.RUM.Parameters(10, "normal")
#' Generate.RUM.Data(Params,m=10,n=5,"normal")
#' Params = Generate.RUM.Parameters(10, "exponential")
#' Generate.RUM.Data(Params,m=10,n=5,"exponential")
Generate.RUM.Data <- function(Params, m, n, distribution) {
  if(distribution == "exponential")
    t(replicate(n, order(rexp(n = m, rate = 1/Params$Mean))))
  else if (distribution == "normal")
    t(replicate(n, order(-rnorm(n = m, mean = Params$Mean, sd = Params$SD))))
  else stop(paste0("Distribution name \"", distribution, "\" not recognized"))
}

Generate.RUM.Multitype.Data <- function(Params, n) {
  K <- Params$K
  m <- Params$m
  getarank <- function() {
    draws <- matrix(rnorm(n = K * m, mean = Params$Mean, sd = Params$SD), nrow = K)
    whichtype <- rmultinom(1, 1, Params$Gamma)
    order(-t(whichtype) %*% draws)
  }
  t(replicate(n, getarank()))
}
################################################################################
# Model Likelihood
################################################################################

#' A faster Likelihood for Plackett-Luce Model
#'
#' @param Data ranking data
#' @param parameter Mean of Exponential Distribution
#' @return log likelihood
#' @export
#' @examples
#' data(Data.Test)
#' parameter = Generate.RUM.Parameters(5, "exponential")
#' Likelihood.PL(Data.Test, parameter)
Likelihood.PL <- function(Data, parameter) {
  gam <- parameter$Mean^-1
  n <- dim(Data)[1]
  m <- dim(Data)[2]
  rank <- Data
  ll <- 0
  
  for(i in 1:n) {
    mj <- sum(rank[i,]>0)
    temp <- gam[rank[i,][rank[i,] > 0]]
    temp <- temp[mj:1]
    ll <- ll+sum(log(temp))-sum(log(cumsum(temp)))
  }
  ll
}

# Helper for Gumbel pdf
pdfgumbel <- function(x,mu) {
  exp(-x+mu)*exp(-exp(-x+mu))
}

#' Likelihood for general Random Utility Models
#'
#' @param Data ranking data
#' @param parameter Mean of Exponential Distribution
#' @param dist exp or norm
#' @param range range
#' @param res res
#' @param race TRUE if data is sub partial, FALSE (default) if not
#' @return log likelihood
#' @export
#' @examples
#' data(Data.Test)
#' parameter = Generate.RUM.Parameters(5, "normal")
#' Likelihood.RUM(Data.Test,parameter, "norm")
Likelihood.RUM <- function(Data, parameter, dist = "exp", range = NA, res = NA, race=FALSE) { 
  if(is.na(range))
    range <- max(abs(parameter$Mean)) + 3 * max(parameter$SD)
  if(is.na(res))
    res <- range / 10000
    
  if(!(dist=="fexp" | dist=="dexp" | dist=="exp" | dist=="norm" | dist=="norm.fixedvariance")) {
    stop(paste0("Distribution name \"", dist, "\" not recognized"))
  }
  
  rank <- Data
  S <- range/res
  if(dist =="fexp" | dist =="dexp" | dist =="exp") {
    x <- (0:S)*res
  }
  else {
    x <- (-S:S)*res
  }

  n <- dim(rank)[1]
  m <- dim(rank)[2]
  
  if(dist=="norm" | dist=="norm.fixedvariance"){
    ll <- 0
    for(i in 1:n){
      if(sum(Data[i,]==0)>0){
        mj <- min(which(rank[i,]==0))-1
        CDF <- matrix(1,1,length(x))
        if(!race){
          for(jt in setdiff(1:m,rank[i,1:mj])){
            CDF <- res*cumsum(dnorm(x,mean=parameter$Mean[jt],sd=parameter$SD[jt]))*CDF
          }
        }
        #mt=min(which(Data[i,]==0))-1
        for(j in mj:1){
          PDF <- dnorm(x,mean=parameter$Mean[rank[i,j]],sd=parameter$SD[rank[i,j]])*CDF
          CDF <- res*cumsum(PDF)
        }
        ll <- log(CDF[length(x)])+ll
        
      }
      if(!(sum(Data[i,]==0)>0)){
        CDF <- matrix(1,1,length(x))
        for(j in m:1){
          PDF <- dnorm(x,mean=parameter$Mean[rank[i,j]],sd=parameter$SD[rank[i,j]])*CDF
          CDF <- res*cumsum(PDF)
        }
        ll <- log(CDF[length(x)])+ll
      }
    }
  }
  
  if(dist=="exp"){
    ll <- 0
    for(i in 1:n){
      if(sum(Data[i,]==0)>0){
        mj <- min(which(rank[i,]==0))-1
        CDF <- matrix(1,1,length(x))
        if(!race){
          for(jt in setdiff(1:m,rank[i,1:mj])){
            CDF <- res*cumsum(dexp(x,rate=parameter$Mean[jt]^-1))*CDF
          }
        }
        
        for(j in mj:1){
          PDF <- dexp(x,rate=parameter$Mean[rank[i,j]]^-1)*CDF
          CDF <- res*cumsum(PDF)
          CDF <- CDF[length(CDF)]-CDF
        }
        ll <- log(CDF[1])+ll
      }
      
      if(!(sum(Data[i,]==0)>0)){
        CDF <- matrix(1,1,length(x))
        for(j in m:1){
          PDF <- dexp(x,rate=parameter$Mean[rank[i,j]]^-1)*CDF
          CDF <- res*cumsum(PDF)
          CDF <- CDF[length(CDF)]-CDF
        }
        ll <- log(CDF[1])+ll
      }
    }
  }
  
  if(dist=="fexp"){
    ll <- 0
    for(i in 1:n){
      if(sum(Data[i,]==0)>0)
      {
        mj <- min(which(rank[i,]==0))-1
        CDF <- matrix(1,1,length(x))
        for(j in 1:mj){
          PDF <- dexp(x,rate=parameter$Mean[rank[i,j]]^-1)*CDF
          CDF <- res*cumsum(PDF)
          CDF <- CDF[length(CDF)]-CDF
        }
        
        if(!race){
          for(jt in setdiff(1:m,rank[i,1:mj]) ){
            PDF <- res*cumsum(dexp(x,rate=parameter$Mean[jt]^-1))*CDF
          }
          CDF <- res*cumsum(PDF)
          CDF <- CDF[length(CDF)]-CDF
        }
        ll <- log(CDF[1])+ll
      }
      
      if(!(sum(Data[i,]==0)>0)){
        CDF <- matrix(1,1,length(x))
        for(j in m:1){
          PDF <- dexp(x,rate=parameter$Mean[rank[i,j]]^-1)*CDF
          CDF <- res*cumsum(PDF)
          CDF <- CDF[length(CDF)]-CDF
        }
        ll <- log(CDF[1])+ll
      }
    }
  }
  
  if(dist=="dexp"){
    x <- (-S:S)*res
    ll <- 0
    for(i in 1:n){
      CDF <- matrix(1,1,length(x))
      for(j in m:1){
        PDF <- pdfgumbel(x,mu=parameter$Mean[rank[i,j]])*CDF
        CDF <- res*cumsum(PDF)
      }
      ll <- log(CDF[length(x)])+ll
    }
  }
  
  ll
}



#' Likelihood for Multitype Random Utility Models
#'
#' @param Data n by m table of rankings
#' @param Estimate Inference object from Estimation function
#' @param dist Distribution of noise (exp or norm)
#' @param race TRUE if data is sub partial, FALSE (default) if not
#' @return log likelihood
#' @export
#' @examples
#' Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
#' Estimate <- Estimation.RUM.MultiType.MLE(Data.Tiny, K=2, iter = 1, dist= "norm")
#' Likelihood.RUM.Multitype(Data.Tiny, Estimate, dist = "norm")
Likelihood.RUM.Multitype <- function(Data, Estimate, dist, race = FALSE) {
  llfull=0
  n <- nrow(Data)
  K <- Estimate$K
  dist <- Estimate$Dist
  for(i in 1:n){
    lli <- 0
    for(k in 1:K){
      if(dist == "exp") {
        ranget <- max(abs(Estimate$Mean[k,]))+3*max(abs(Estimate$Mean[k,]))^.5
        ll <- Likelihood.RUM(matrix(as.numeric(Data[i,]), nrow=1),parameter = list(Mean=Estimate$Mean[k,]),dist = "fexp", range = ranget, res = ranget * .0001, race = race)
      }
      if(dist == "norm" | dist == "norm.fixedvariance"){
        ranget <- max(abs(Estimate$Mean[k,]))+3*max(Estimate$SD[k,])
        ll <- Likelihood.RUM(matrix(as.numeric(Data[i,]), nrow=1),parameter=list(Mean=Estimate$Mean[k,],SD=Estimate$SD[k,]),dist,range=ranget,res=ranget*.0001, race = race)
      }
      lli <- exp(ll)*Estimate$Gamma[k]+lli
    }
    llfull <- llfull+log(lli)
  }
  llfull
}


################################################################################
## Estimation Methods MLE: MCEM( MCEM.Ag() )
################################################################################

#' Performs parameter estimation for the Plackett-Luce model using an Minorize Maximize algorithm
#' 
#' @param Data data in either partial or full rankings (Partial rank case works for settings like car racing)
#' @param iter number of MM iterations to run
#' @return list of estimated means (Gamma) and the log likelihoods
#' @export
#' @examples
#' data(Data.Test)
#' Estimation.PL.MLE(Data.Test)
Estimation.PL.MLE <- function(Data, iter = 10) {
  rank <- Data
  t0 <- proc.time()
  m <- dim(rank)[2]
  n <- dim(rank)[1]  
  
  GamaTotal <- matrix(0,iter,m)
  LTotal <- matrix(0,1,iter)
  
  M <- matrix(0,1,n)
  for(i in 1:n){
    if(sum(rank[i,]==0)){
      M[i] <- min(which(rank[i,]==0))-1
    }
    if(!sum(rank[i,]==0)){
      M[i] <- m
    }
  }
  
  W <- matrix(0,1,m)
  #nominator
  for(t in 1:m){
    W[t] <- 0
    for(j in 1:n){
      W[t] <- (rank[j,M[j]]==t)+W[t]
    }
    W[t] <- n-W[t]
  }
  
  gam <- matrix(1,1,m)
  for(iteration in 1:iter){
    gamtemp <- gam 
    
    for(t in 1:m){     
      denom <- 0.1
      for(j in 1:n){
        for(i in 1:(M[j]-1)){
          delta <- 0
          delta <- sum(rank[j,i:M[j]]==t)
          denomt3 <- 0
          for( s in i:M[j]){
            denomt3 <- denomt3+gam[rank[j,s]]
          }
          denom <- delta/denomt3+denom   
        }
      }
      gamtemp[t] <- W[t]/denom    
    }
    
    gam=gamtemp
    GamaTotal[iteration,]=gam/sum(gam)
    
    ll=0
    for(i in 1:n){
      mj <- sum(rank[i,]>0)
      temp <- gam[rank[i,][rank[i,] > 0]]
      temp <- temp[mj:1]
      ll <- ll+sum(log(temp))-sum(log(cumsum(temp)))
    }
    LTotal[iteration] <- ll
    print(paste0("Finished ", iteration, "/", iter))
  } 
  
  means <- (GamaTotal[iter,]^-1) /sum(GamaTotal[iter,]^-1)
  
  t <- proc.time() - t0
   
  list(m = m,
       order = order(means),
       Mean = means,
       SD = means,
       LL = LTotal, 
       Time = t,
       AverageLogLikelihood=t(LTotal[iter]/n),
       Parameters = convert.vector.to.list(means))
}

################################################################################
## Estimation Methods MCEM( MCEM.Ag() )
################################################################################

#' Performs parameter estimation for a Random Utility Model with different noise distributions
#' 
#' This function supports RUMs 
#' 1) Normal
#' 2) Normal with fixed variance (fixed at 1)
#' 3) Exponential (top k setting like Election)
#' 
#' @param Data data in either partial or full rankings
#' @param iter number of EM iterations to run
#' @param dist underlying distribution. Can be "norm", "norm.fixedvariance", "exp"
#' @param race indicator that each agent chose a random subset of alternatives to compare
#' @return parameters of the latent RUM distributions
#' @export
#' @examples
#' Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
#' Estimation.RUM.MLE(Data.Tiny, iter = 2, dist="norm")
Estimation.RUM.MLE = function(Data, iter = 10, dist, race = FALSE)
{
  if(!(dist=="dexp" | dist=="exp" | dist=="norm" | dist == "norm.fixedvariance" ))
    stop(paste0("Distribution name \"", dist, "\" not recognized"))
  
  # calculating the dimensions
  dims <- dim(Data)
  n <- dims[1]
  m <- dims[2]
  
  #flipping the data for exponential distribution
  if(dist=="exp")
  {
    DataL <- Data
    Data <- Data[, m:1]
  }
  
  t0 <- proc.time()
  
  #initialization of mean and variance
  Delta <- exp(rnorm(m))
  Variance <- exp(rnorm(m))
  if(dist=="norm" || dist == "norm.fixedvariance"){
    parameter <- list(Mu=Delta, Var=Variance)
    dist2 <- "norm"
  }
  if(dist=="exp"){
    parameter <- list(Mu=Delta)
    dist2 <- "exp"
  }
  
  ##############
  MM=matrix(0,iter,m)
  VV=matrix(0,iter,m)
  ll=matrix(0,1,iter)
  
  #EM Algorithm
  ##############
  ##############
  for(j in 1:iter)
  {
    S=2000+300*j
    print(paste0("Iteration ", j,"/",iter))
    
    U=matrix(0,n,m)
    U2=matrix(0,n,m)
    
    #E-Step-Parallelizable:
    ##############
    ##############
    for(k in 1:n)
    {
      initial=matrix(0,1,m)
      print(k)
      initial[Data[k,][Data[k,]>0]]=sort(runif(sum(Data[k,]>0)),decreasing=TRUE)
      Temp=GibbsSampler(S,initial,pi=Data[k,],dist2,parameter, race)
      U[k,]=Temp$M1
      if(race){U[k,which(Data[k,]==0)]=Delta[which(Data[k,]==0)]}
      U2[k,]=Temp$M2
    }
    
    # Mstep
    ##############
    ##############
    
    Delta=1/n*colSums(U)
    if(dist != "exp"){
      Delta=Delta-Delta[1]+1
    }
    if(dist == "norm.fixedvariance") {
      Variance <- matrix(.1,1,m)
    }
    if(dist == "exp") {
      Variance <- Delta^2
    }
    if(dist == "norm"){
      Variance <- (1/n*colSums(U2-2* U * (matrix(1,n,1) %*% parameter$Mu)+(matrix(1,n,1) %*% parameter$Mu)^2))
      #Variance <- Variance / sum(Variance)
      #Variance=(1/n*colSums(U2)-Delta^2)
      #Variance[1]=0.1
    }
        
    if(dist=="norm" || dist == "norm.fixedvariance"){
      parameter <- list(Mu=Delta,Var=Variance)
    }
    if(dist=="exp"){
      parameter <- list(Mu=Delta/sum(Delta))
    }
    
    MM[j,]=Delta
    VV[j,]=Variance
    print(parameter$Mu)
    print(Variance)   
  }  
 
  RT <- sort(MM[iter,],decreasing=TRUE,index=TRUE)
  t <- proc.time() - t0
  
  params <- rep(list(list()), m)
  for(i in 1:m) {
    if(dist == "exp") params[[i]]$Mean <- 1/Delta[i]
    else params[[i]]$Mean <- Delta[i]
    params[[i]]$SD <- Variance[i]^.5
  }
  ###
  
  
  if(dist == "exp") ordering <- order(MM[iter,])
  else ordering <- order(-MM[iter,])
  list(m = m,
       order = ordering,
       Aggregated.Rank=RT$ix,
       Mean=Delta,
       SD=Variance^.5,
       Time = t,
       Parameters = params) 
  
}

#' Performs parameter estimation for a Multitype Random Utility Model
#' 
#' This function supports RUMs 
#' 1) Normal
#' 2) Normal with fixed variance (fixed at 1)
#' 3) Exponential
#' 
#' @param Data data in either partial or full rankings
#' @param K number of components in mixture distribution
#' @param iter number of EM iterations to run
#' @param dist underlying distribution. Can be "norm", "norm.fixedvariance", "exp"
#' @param ratio parameter in the algorithm that controls the difference of the starting points, the bigger the ratio the more the distance
#' @param race TRUE if data is sub partial, FALSE (default) if not
#' @return results from the inference
#' @export
#' @examples
#' Data.Tiny <- matrix(c(1, 2, 3, 3, 2, 1, 1, 2, 3), ncol = 3, byrow = TRUE)
#' Estimation.RUM.MultiType.MLE(Data.Tiny, K=2, iter = 3, dist= "norm.fixedvariance")
Estimation.RUM.MultiType.MLE = function(Data, K = 2, iter = 10, dist, ratio = .2, race=FALSE)
{
  if(!(dist=="dexp" | dist=="exp" | dist=="norm" | dist == "norm.fixedvariance" ))
    stop(paste0("Distribution name \"", dist, "\" not recognized"))
  
  # calculating the dimensions
  dims <- dim(Data)
  n <- dims[1]
  m <- dims[2]
  
  #flipping the data for exponential distribution
  if(dist=="exp")
  {
    DataL <- Data
    Data <- Data[, m:1]
  }
  
  t0 <- proc.time()
  
  #initialization of mean and mixture probabilities
  initiateMean <- matrix(exp(rnorm(m*K)),K,m)
  for(k in 1:K) {
    initiateMean[k,] <-  initiateMean[1,] + runif(m)*ratio
    initiateMean[k,1] <- 1
  }
  
  #initiateGamma <- runif(K)
  #initiateGamma <- initiateGamma/sum(initiateGamma)
  initiateGamma <- matrix(1/K,1,K)
  
  initiateZ <- matrix(sample(1:K, size=n, replace = TRUE, prob = NULL),1,n)
  
  Delta <- initiateMean
  if(dist!="exp"){
    Delta[,1] <- matrix(1,1,K)
  }
  Variance <- matrix(.1,K,m)
  Z <-initiateZ
  gamma=initiateGamma
  
  ##############
  if(dist=="norm" | dist == "norm.fixedvariance"){
    parameter <- list(Mu=Delta, Var=Variance)
    dist2 <- "norm"
  }
  if(dist=="exp"){
    parameter <- list(Mu=Delta, Var=Delta^2)
    dist2 <- "exp"
  }
  
  #EM Algorithm
  ############## 
  ############## 
  
  for(j in 1:iter){
    #setting number of samples to increase by iterations
    S <- 200+300*j
    #S <- 1000+100*j
    #S = 100
    print(sprintf("Iteration %d",j))
  
    z <- matrix(0,1,n)
    ESTEP1 <- matrix(0,K,m)
    ESTEP2 <- matrix(1,n,K)
    
    #E-Step-Parallelizable:
    #############
    #############
    
    initial <- matrix(0,1,m)
    #print(Delta)
    #print(gamma)
    
    for(i in 1:n){
      initial[Data[i,][Data[i,]>0]] <- sort(runif(sum(Data[i,]>0)),decreasing=TRUE)
      X <- initial
      for(t1 in 1:S) { 
        ZProb <- gamma*PdfModel(X,parameter,dist)
        ZProb <- ZProb/sum(ZProb)
        zi <- sample(1:K, size = 1, replace = TRUE, prob = ZProb)
        XTemp <- GibbsSampler(3,X,pi=Data[i,],dist2,parameter=list(Mu=Delta[zi,],Var=Variance[zi,]),race)
        X <- XTemp$M1
        #X <- XTemp$X
        #print(X)
        ESTEP1[zi,] <- ESTEP1[zi,] + X 
        ESTEP2[i,zi] <- ESTEP2[i,zi] + 1
      }
    }
    
    #print(ESTEP2)
    # Mstep
    ##############
    ##############
    
    for(k in 1:K){
      Delta[k,] <- 1/(sum(ESTEP2[,k]))*ESTEP1[k,]
    }
    
    if(dist!="exp"){
      Delta[,1] <- matrix(1,1,K)
      #Delta <- Delta-Delta[1,1]+1
      #Variance=matrix(.1,K,m)
    }
    
    for(k in 1:K){
      gamma[k]=sum(ESTEP2[,k])
      #gamma[k] <- sum(ESTEP2[,k]+1)
    }
    gamma <- gamma/sum(gamma)
    
    #gamma=gamma+matrix(1/K,1,K)
    #gamma = gamma/sum(gamma)
    
    if(dist=="norm" | dist == "norm.fixedvariance") parameter <- list(Mu=Delta,Var=Variance)
    if(dist=="exp") { 
      for(k in 1:K) Delta[k,] <- Delta[k,]/sum(Delta[k,])
      parameter <- list(Mu=Delta)
    }
  }
    
    #############
    #############
 
  t <- proc.time()-t0
  
  params <- rep(list(list()), m)
  for(i in 1:m) {
    params[[i]]$Mean <- Delta[,i]
    params[[i]]$SD <- Variance[,i]^.5
    params[[i]]$Gamma <- gamma
  }
  
  ###Computing loglikelihood
  #LL = Likelihood.RUM.Multitype(DataL, Estimate =  list(Mean = Delta, SD = Variance[,i]^.5, Gamma = gamma, K = K, Dist = dist), race = race)
  
  if(dist == "exp") ordering <- order(gamma %*% Delta)
  else ordering <- order(-gamma %*% Delta)
    
  list(m = m, order = ordering, Dist = dist, K = K, Mean=Delta, SD=Variance^.5, Gamma=gamma, Personal=ESTEP2, 
       Time=t, Parameters = params) 
  
}



# The PDF of a set of independent normal random variables
# 
# 
# @param X vector of scalar values
# @param parameter list containing vetor of means for normal distributions, variances are set to one
# @return value of the PDF for vector X given the parameter
# @export
PdfModel <- function(X, parameter, dist = "norm") {
  m <- length(X)
  K <- dim(parameter$Mu)[1]
  out <- matrix(0,1,K)
  
  if(dist == "norm" || dist == "norm.fixedvariance"){
    for(k in 1:K)
      out[k] <- prod(dnorm(X,mean=parameter$Mu[k,],sd=parameter$Var[k]^.5))
  }
  if(dist == "exp"){
    for(k in 1:K)
      out[k] <- prod(dexp(X,rate=parameter$Mu[k,]^-1))
  }
  out
}


#' Performs parameter estimation for a Generalized Random Utility Model with user and alternative characteristics
#' 
#' This function supports RUMs 
#' 1) Normal with fixed variance (fixed at 1)
#' 
#' 
#' @param Data data in either partial or full rankings
#' @param X user characteristics
#' @param Z alternative characteristics
#' @param iter number of iterations to run algorithm
#' @param dist choice of distribution
#' @param din initialization of delta vector
#' @param Bin intialization of B matrix
#' @return results from the inference
#' @export
#' @examples
#' #data(Data.Test)
#' #Data.X= matrix( runif(15),5,3)
#' #Data.Z= matrix(runif(10),2,5)
#' #Estimation.GRUM.MLE(Data.Test, Data.X, Data.Z, iter = 3, dist = "norm", 
#' #din=runif(5), Bin=matrix(runif(6),3,2))
Estimation.GRUM.MLE <- function(Data, X, Z, iter, dist, din, Bin) {
  T0=proc.time() 
  n <- dim(Data)[1]
  m <- dim(Data)[2]
  
  U=X%*%Bin%*%Z+t(matrix(din,m,1)%*%matrix(1,1,n))
  #U2=matrix(abs(rnorm(n*m,0,1)),n,m)
  U2=matrix(1,n,m)
  
  for(j in 1:iter)
  {
    S=1000+500*j
    print(j)
    ##########################E-Step-Parallelizable#########################
    for(k in 1:n)
    {
      initial=matrix(0,1,m)
      initial[Data[k,which(Data[k,]>0)]]=sort(runif(sum(Data[k,]>0)),decreasing=TRUE)
      Temp=GibbsSampler(S,initial,pi=Data[k,],dist,parameter=list(Mu=U[k,],Variance=U2[k,]))
      U[k,]=Temp$M1
      #U2[k,]=Temp$M2
    }
    ############################### Mstep###################################
    
    psi=RegUXZ(U,X,Z)
    deltah=psi$deltah
    Bh=psi$Bh
    
    deltah[1]=0
    U=X%*%Bh%*%Z+t(matrix(deltah,m,1)%*%matrix(1,1,n))
    #U2[,1]=1
    ########################################################################
    
  }  
  
  DT=proc.time()-T0
  
  list(m = m, Par=psi,Time=DT) 
  
}

RegUXZ <- function(U,X,Z) {
  n <- dim(U)[1]
  m <- dim(U)[2]
  K <- dim(X)[2]
  L <- dim(Z)[1]
  
  Uv <- matrix(U,m*n,1)
  Xv <- matrix(0,n*m,K*L+m)
  
  for(row in 1:(n*m)) {
    j <- ceiling(row/n)
    i <- row-(j-1)*n
    Xv[row,1:(K*L)] <- matrix(X[i,]%*%t(Z[,j]),1,K*L)
    Xv[row,K*L+j] <- 1
  }
  
  out <- glm.fit(Xv,Uv,intercept = FALSE)
  Bh <- matrix(out$coefficients[1:(K*L)],K,L)
  deltah <- out$coefficients[(K*L+1):(K*L+m)]
  list(Bh=Bh,deltah=deltah,MU=X%*%Bh%*%Z+t(matrix(deltah,m,1)%*%matrix(1,1,n)))
}


################################################################################
## GMM
################################################################################

#' GMM Method for estimating Plackett-Luce model parameters
#' 
#' @param Data.pairs data broken up into pairs
#' @param m number of alternatives
#' @param prior magnitude of fake observations input into the model
#' @param weighted if this is true, then the third column of Data.pairs is used as a weight for that data point
#' @return Estimated mean parameters for distribution of underlying exponential
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' Estimation.PL.GMM(Data.Test.pairs, 5)
Estimation.PL.GMM <- function(Data.pairs, m, prior = 0, weighted = FALSE) { 
  
  t0 <- proc.time()
  
  transition.matrix <- generateC(Data.pairs, m, weighted = weighted, prior, normalized = FALSE)
  
  diag(transition.matrix) <- diag(transition.matrix) - min(diag(transition.matrix))
  transition.matrix <- transition.matrix / rowSums(transition.matrix)
  
  first.eigenvector <- Re(eigen(t(transition.matrix))$vectors[,1])
  stationary.probability <- first.eigenvector / sum(first.eigenvector)
  means <- stationary.probability / sum(stationary.probability)
  
  # this is the estimated means
  list(m = m,
       order = order(means),
       Mean = means,
       SD = means,
       Time = proc.time() - t0,
       Parameters = convert.vector.to.list(stationary.probability))
}

#' GMM Method for Estimating Random Utility Model wih Normal dsitributions
#' 
#' @param Data.pairs data broken up into pairs
#' @param m number of alternatives
#' @param iter number of iterations to run
#' @param Var indicator for difference variance (default is FALSE)
#' @param prior magnitude of fake observations input into the model
#' @return Estimated mean parameters for distribution of underlying normal (variance is fixed at 1)
#' @export
#' @examples
#' data(Data.Test)
#' Data.Test.pairs <- Breaking(Data.Test, "full")
#' Estimation.Normal.GMM(Data.Test.pairs, 5)
Estimation.Normal.GMM <- function(Data.pairs, m, iter=1000, Var = FALSE, prior = 0) {
    
  t0 <- proc.time()
  
  sdhat <- matrix(1,1,m)
  muhat <- matrix(1,1,m)
  C <- generateC(Data.pairs, m, prior)

  if(!Var) {
    for(iter in 1:iter) {
      alpha <- iter^-1
      muhat <- muhat + alpha *rowSums(exp(-delta(muhat)^2/4)*(C - f(muhat)))
      muhat <- muhat - min(muhat)
    }
    #print("Error = ")
    print(sum((C - f(muhat))^2)^.5)
  }
  
  if(Var) {
    for(iter in 1:iter) {
      alpha <- iter^(-1)
      muhat <- muhat + alpha *rowSums(exp(-delta(muhat,sdhat,Var=TRUE)^2/4)*(C - f(muhat,sdhat,Var=TRUE)))
      sdhat <- abs(sdhat - alpha *rowSums(VarMatrix(sdhat)^(-2)*exp(-delta(muhat,sdhat,Var=TRUE)^2/4)*(C - f(muhat,sdhat,Var=TRUE))))
      muhat <- muhat - min(muhat)
      sdhat[1] <- 1
    }
    #print("Error = ")
    print(sum((C - f(muhat,sdhat,Var=TRUE))^2)^.5)
  }
  
  t <- proc.time() - t0
  
  params <- rep(list(list()), m)
  for(i in 1:m) {
    params[[i]]$Mean <- muhat[1, i]
    params[[i]]$SD <- sdhat[1, i]
  }
  
  list(m = m, order = order(-muhat[1,]), Mean = muhat[1,], SD = sdhat[1,], Time = t, Parameters = params)
}