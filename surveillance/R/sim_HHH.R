###################################################
### chunk number 1: 
###################################################
############################################
# Simulates multivariate count data based on the model described in Held et.al (2005)
# Note: trend is omitted
######################################
simHHH.default <- function(model=NULL,control=list(coefs=list(alpha=1, gamma=0, delta=0, lambda=0,
                                         phi=NULL,psi=NULL,period=52),
                           neighbourhood=NULL,population=NULL,start=NULL),
                           length){

   #################################################
   #Help functions
   ################################################

   # draws n random numbers from a NB(mean, psi) distribution
   rNB<-function(n,mean,size=control$coefs$psi){
      rnbinom(n, mu=mean, size=size)
   }

   # returns formula for the seasonal part of \nu_t
   formulaSeason <- function(mod="~-1",S=1,period){
    for(i in 1:S){
      mod <- paste(mod,"+sin(",2*i,"*pi*t/",period,")+cos(",2*i,"*pi*t/",period,")",sep="")
    }
    return(as.formula(mod))
  }

  # sum of all neighbours
  # params: x - vector with counts
  #         nhood - adjacency matrix, 0= no neighbour
  # returns a vector with the sum of "neighbouring counts" for all areas
  sumN <- function (x,nhood) {
    n<- length(x)
    if(any(nhood>0)){
      nhood <- nhood >0
      res <- sapply(1:n,function(i) sum(x[nhood[i,]]))
    } else {
      res<- rep(0,n)
    }
    return(res)
  }

  ##################################################

   ##################################
   # set default values if not specified
   #####################################

   if(is.null(control$coefs$alpha))
     stop("alpha needs to be specified")
   nAreas <- length(control$coefs$alpha)

   # define neighbourhood-matrix, assume there are no neighbours
   if(is.null(control$neighbourhood))
     control$neighbourhood <- matrix(0,nAreas,nAreas)

   # set population (i.e. n_i,t) to 1 if not specified
   if(is.null(control$population)){
     control$population <- matrix(1, ncol=nAreas,nrow=length)
   } else {
     #assumption: n_i,t = n_i
     pop <-control$population[1,]
     control$population <- matrix(pop,ncol=nAreas,nrow=length,byrow=TRUE)
   }

   #determine number of seasons
   if(is.null(control$coefs$gamma)){
     control$coefs$gamma <-0
     control$coefs$delta <- 0
     S <- 1
   } else {
     if(length(control$coefs$gamma) != length(control$coefs$delta))
       stop("gamma and delta must have the same length")
     S <- length(control$coefs$gamma)
   }
   if(is.null(control$coefs$period))
     control$coefs$period <- 52

   # is there a autoregressive (epidemic) part
   if(is.null(control$coefs$lambda)){
     control$coefs$lambda <- 0
   }
   if(is.null(control$coefs$phi)){
     control$coefs$phi <- 0
   }
   if(!is.null(control$start)){
     if(length(control$start)!=nAreas)
       stop("wrong dimension of start\n")
   }

   # simulate from Poisson or NegBin model
   if(is.null(control$coefs$psi)){
     rdistr<-rpois
   } else{
      rdistr<-rNB
   }

   # computation of seasonal part of nu_i,t:
    season <- model.frame(formula=formulaSeason(S=S,period=control$coefs$period),
                           data=data.frame("t"=1:length))
    #rearrange the sinus and cosinus parts
    season <- season[,c(seq(1,2*S,by=2),seq(2,2*S,by=2))]
    # this computes \sum_{s=1}^S [gamma_s*sin(omega_s*t) + delta_s*cos(omega_s*t) ]
    season<- as.matrix(season)%*%c(control$coefs$gamma,control$coefs$delta)

   # compute endemic part: nu_t =  exp( alpha_i + season_t )
    nu<-exp(sapply(1:nAreas,function(i) control$coefs$alpha[i]+season))

    # initialize matrices for the mean mu_i,t and the simulated data x_i,t
    # x_i,0 is set to the mean of n_it*\nu_it
    mu <- matrix(0,ncol=nAreas,nrow=length)
    x <- matrix(0,ncol=nAreas,nrow=length+1)

    x[1,] <- ifelse(is.null(control$start),colMeans(control$population*nu),start)
    #print(x[1,])

    if(control$coefs$lambda == 0 && control$coefs$phi ==0){
      mu <- control$population*nu
      x <- matrix(rdistr(nAreas*(length+1),mu),ncol=nAreas,byrow=FALSE)
    } else {
      # simulate data
      for(t in 1:length){
        #mu_i,t = lambda*x_i,t-1 +phi*\sum_j~i x_j,t-1
        mu[t,] <- control$coefs$lambda *x[t,] + control$coefs$phi*sumN(x[t,],control$neighbourhood) + control$population[t,]*nu[t,]
        x[t+1,] <- rdistr(nAreas,mu[t,])
      }
    }
    #remove first time point
    dp <- create.disProg(week=1:length,observed=x[-1,],state=rep(0,length),
                         neighbourhood=control$neighbourhood, populationFrac=control$population)
    return(list(data=dp,mean=mu,endemic=control$population*nu,coefs=control$coefs))
}

##################
simHHH <- function(model,control,length){
  UseMethod("simHHH")
}

################################
# simulates data using the estimated parameter values of a model fitted with algo.hhh
# Note: NO trend
simHHH.ah <- function(model,control=model$control, length){
  #hoehle: removed this to make simHHH.ah consistent with simHHH.default
  # control <- model$control

  #number of areas
  nAreas <- ncol(model$disProgObj$observed)
  #number of seasons
  S <- control$nseason

  cntrl <- list(lambda=NULL,phi=NULL,gamma=NULL,delta=NULL,
                 psi=NULL,period=model$control$period)

 #extract coefficients
 coefs <- coef(model)
  if(control$neighbours)
    cntrl$phi <- coefs["phi"]
  if(control$negbin)
    cntrl$psi <- coefs["psi"]
  if(control$lambda)
    cntrl$lambda <- coefs["lambda"]

  if(S > 0){
    cntrl$gamma <- coefs[paste("gamma",1:S,sep="")]
    cntrl$delta <- coefs[paste("delta",1:S,sep="")]
  }

  cntrl$alpha <- coefs[paste("alpha",1:nAreas,sep="")]

   result <- simHHH(length,control=list(coefs=cntrl,
                                  neighbourhood=model$disProgObj$neighbourhood,
                                  populationFrac=model$disProgObj$populationFrac
                                  ))
   return(result)
}


