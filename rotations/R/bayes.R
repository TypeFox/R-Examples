#' MCMC for rotation data
#'
#' Use non-informative Bayesian methods to infer about the central orientation and concentration 
#' parameter for a sample of rotations.
#'
#' The procedures detailed in \cite{bingham2009b} and \cite{bingham2010} are implemented to obtain
#' draws from the posterior distribution for the central orientation and concentration parameters for 
#' a sample of 3D rotations.  A uniform prior on SO(3) is used for the central orientation and the
#' Jeffreys prior determined by \code{type} is used for the concentration parameter.  
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type Angular distribution assumed on R.  Options are \code{\link{Cayley}}, \code{\link{Fisher}} or \code{\link{Mises}}
#' @param S0 initial estimate of central orientation
#' @param kappa0 initial estimate of concentration parameter
#' @param tuneS central orientation tuning parameter, concentration of proposal distribution
#' @param tuneK concentration tuning parameter, standard deviation of proposal distribution
#' @param burn_in number of draws to use as burn-in
#' @param m number of draws to keep from posterior distribution
#' @return  list of 
#' \itemize{
#'          \item \code{S} Draws from the posterior distribution for central orientation S
#'          \item \code{kappa} Draws from the posterior distribution for concentration parameter kappa
#'          \item \code{Saccept} Acceptance rate for central orientation draws
#'          \item \code{Kaccept} Acceptance rate for concentration draws
#'          }
#' @@cite bingham2009b bingham2010
#' @export
#' @examples
#' #Not run due to time constraints
#' \dontrun{
#' Rs <- ruars(20, rfisher, kappa = 10)
#' draws <- MCMCSO3(Rs, type = "Fisher", S0 = mean(Rs), kappa0 = 10, tuneS = 5000, 
#'                  tuneK = 1,burn_in = 1000, m = 5000)}

MCMCSO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  UseMethod("MCMCSO3")
}


#' @rdname MCMCSO3
#' @method MCMCSO3 SO3
#' @export 

MCMCSO3.SO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  if(type %in% c("Cayley","cayley")){
    
    lpangle <- 1
    
  }else if(type %in% c("Fisher","fisher")){
    
    lpangle <- 2
    
  }else if(type %in% c("Mises","mises")){
    
    lpangle <- 3
    
  }else{
    stop("Invalid choise of type: please choose Cayley, Fisher or Mises.")
  }
  
  listRes<-both_MCMC_CPP(x,S0, kappa0,tuneS,tuneK,burn_in,m, lpangle)
  class(listRes$S)<-"SO3"
  #listRes$S<-as.SO3(listRes$S)
  
  return(listRes)
}

#' @rdname MCMCSO3
#' @method MCMCSO3 Q4
#' @export 

MCMCSO3.Q4<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  Rs<-as.SO3(x)
  S0<-as.SO3(matrix(S0,3,3))
  SO3Res<-MCMCSO3(Rs,type,S0,kappa0,tuneS,tuneK,burn_in,m)
  Q4Res<-list(Q=as.Q4(SO3Res$S),kappa=SO3Res$kappa,Qaccept=SO3Res$Saccept,Kaccept=SO3Res$Kaccept)
  return(Q4Res)
  
}

#' Bayes credible regions
#'
#' Find the radius of a \eqn{100(1-\alpha)}\% credible region for the central orientation and concentration parameter using 
#' non-informative Bayesian methods.
#'
#' Compute the radius of a \eqn{100(1-\alpha)}\% credible region for the central orientation and concentration parameter
#' as described in \cite{bingham2009b} and \cite{bingham2010}.  The posterior mode is returned along with the radius
#' of the credible region centered at the posterior mode.
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type Angular distribution assumed on R.  Options are \code{\link{Cayley}}, \code{\link{Fisher}} or \code{\link{Mises}}
#' @param S0 initial estimate of central orientation
#' @param kappa0 initial estimate of concentration parameter
#' @param tuneS central orientation tuning parameter, concentration of proposal distribution
#' @param tuneK concentration tuning parameter, standard deviation of proposal distribution
#' @param burn_in number of draws to use as burn-in
#' @param m number of draws to keep from posterior distribution
#' @param alp alpha level desired, e.g. 0.05 or 0.10.
#' @return  list of 
#'          \itemize{
#'          \item \code{Shat,Qhat} Mode of the posterior distribution for the central orientation S
#'          \item \code{Radius} Radius of the credible region centered at the posterior mode
#'          }
#' @seealso \code{\link{fisheretal}}, \code{\link{prentice}}, \code{\link{chang}}, \code{\link{zhang}}
#' @@cite bingham2009b bingham2010
#' @export
#' @examples
#' #Not run due to time constraints
#' \dontrun{
#' Rs <- ruars(20, rvmises, kappa = 10)
#' 
#' #Compare the region size of the moment based theory mean estimator to the 
#' #Bayes region.
#' 
#' region(Rs, method = "direct", type = "theory", estimator = "mean", alp=0.1, m = 100)
#' bayesCR <- region(Rs, type = "Mises", method = "Bayes", estimator = "mean", S0 = mean(Rs),
#'                    kappa0 = 10, tuneS = 5000, tuneK = 1, burn_in = 1000, alp = .01, m = 5000)
#'                    
#' bayesCR$Radius       #Region size is give by "Radius"
#' bayesCR$Shat         #The Bayes region is centered around the posterior mode: "Shat"}


bayesCR<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000,alp=0.1){
  UseMethod("bayesCR")
}


#' @rdname bayesCR
#' @method bayesCR SO3
#' @export 

bayesCR.SO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000,alp=0.1){
  
  if(type %in% c("Cayley","cayley")){
    
    #lpangle <- lpcayley
    lpangle <- 1
    
  }else if(type %in% c("Fisher","fisher")){
    
    #lpangle <- lpfisher
    lpangle <- 2
    
  }else if(type %in% c("Mises","mises")){
    
    #lpangle <- lpvmises
    lpangle <- 3
    
  }else{
    stop("Invalid choise of type: please choose Cayley, Fisher or Mises.")
  }
  
  listRes<-both_MCMC_CPP(x,S0, kappa0,tuneS,tuneK,burn_in,m, lpangle)
  Sdraws<-as.SO3(listRes$S)
  Shat<-mean(Sdraws)
  rs<-rot.dist(Sdraws,Shat)
  
  return(list(Shat=Shat,Radius=quantile(rs,1-alp)))
}

#' @rdname bayesCR
#' @method bayesCR Q4
#' @export 

bayesCR.Q4<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000,alp=0.1){
  
  Rs<-as.SO3(x)
  S0<-as.SO3(matrix(S0,3,3))
  SO3Res<-bayesCR(Rs,type,S0,kappa0,tuneS,tuneK,burn_in,m,alp)
  Q4Res<-list(Qhat=as.Q4(SO3Res$Shat),Radius=SO3Res$Radius)
  return(Q4Res)
  
}

#' Parameter estimates based on non-informative Bayes
#'
#' Use non-informative Bayes to estimate the central orientation and concentration parameter of a sample of rotations.
#'
#' The procedures detailed in \cite{bingham2009b} and \cite{bingham2010} are implemented to obtain
#' draws from the posterior distribution for the central orientation and concentration parameters for 
#' a sample of 3D rotations.  A uniform prior on SO(3) is used for the central orientation and the
#' Jeffreys prior determined by \code{type} is used for the concentration parameter.  
#'
#' @param x \eqn{n\times p}{n-by-p} matrix where each row corresponds to a random rotation in matrix (\eqn{p=9}) or quaternion (\eqn{p=4}) form.
#' @param type Angular distribution assumed on R.  Options are \code{\link{Cayley}}, \code{\link{Fisher}} or \code{\link{Mises}}
#' @param S0 initial estimate of central orientation
#' @param kappa0 initial estimate of concentration parameter
#' @param tuneS central orientation tuning parameter, concentration of proposal distribution
#' @param tuneK concentration tuning parameter, standard deviation of proposal distribution
#' @param burn_in number of draws to use as burn-in
#' @param m number of draws to keep from posterior distribution
#' @return  list of 
#'          \itemize{
#'          \item \code{Shat} Mode of the posterior distribution for the central orientation S
#'          \item \code{kappa} Mean of the posterior distribution for the concentration kappa
#'          }
#' @seealso \code{\link{mean.SO3}}, \code{\link{median.SO3}}
#' @@cite bingham2009b bingham2010
#' @export
#' @examples
#' Rs <- ruars(20, rvmises, kappa = 10)
#' 
#' Shat <- mean(Rs)               #Estimate the central orientation using the projected mean
#' rotdist.sum(Rs, Shat, p = 2)   #The projected mean minimizes the sum of squared Euclidean
#' rot.dist(Shat)                 #distances, compute the minimized sum and estimator bias 
#' 
#' #Estimate the central orientation using the posterior mode (not run due to time constraints) 
#' #Compare it to the projected mean in terms of the squared Euclidean distance and bias
#' \dontrun{
#' ests <- bayes.mean(Rs, type = "Mises", S0 = mean(Rs), kappa0 = 10, tuneS = 5000,
#'                    tuneK = 1, burn_in = 1000, m = 5000)
#'                    
#' Shat2 <- ests$Shat             #The posterior mode is the 'Shat' object
#' rotdist.sum(Rs, Shat2, p = 2)  #Compute sum of squared Euclidean distances
#' rot.dist(Shat2)                #Bayes estimator bias}

bayes.mean<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  UseMethod("bayes.mean")
}


#' @rdname bayes.mean
#' @method bayes.mean SO3
#' @export 

bayes.mean.SO3<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  if(type %in% c("Cayley","cayley")){
    
    #lpangle <- lpcayley
    lpangle <- 1
    
  }else if(type %in% c("Fisher","fisher")){
    
    #lpangle <- lpfisher
    lpangle <- 2
    
  }else if(type %in% c("Mises","mises")){
    
    #lpangle <- lpvmises
    lpangle <- 3
    
  }else{
    stop("Invalid choise of type: please choose Cayley, Fisher or Mises.")
  }
  
  listRes<-both_MCMC_CPP(x,S0, kappa0,tuneS,tuneK,burn_in,m, lpangle)
  postRes<-list(Shat=mean.SO3(listRes$S),kappa=mean(listRes$kappa))
  
  return(postRes)
}

#' @rdname bayes.mean
#' @method bayes.mean Q4
#' @export 

bayes.mean.Q4<-function(x,type,S0,kappa0,tuneS,tuneK,burn_in,m=5000){
  
  Rs<-as.SO3(x)
  S0<-as.SO3(matrix(S0,3,3))
  SO3Res<-MCMCSO3(Rs,type,S0,kappa0,tuneS,tuneK,burn_in,m)
  Q4Res<-list(Qhat=as.Q4(mean(SO3Res$S)),kappa=mean(SO3Res$kappa))
  return(Q4Res)
  
}
