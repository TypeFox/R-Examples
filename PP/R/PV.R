#' Draw Plausible values
#' 
#' This function draws \code{npv} plausible values for each person from their posterior density.
#' 
#' 
#' @param estobj An object which originates from using \code{PP_4pl()}, \code{PP_gpcm()} or \code{PPall()}. EAP estimation is strongly recommended (\code{type = "eap"}), when plausible values are drawn afterwards, because the EAP estimate is used as starting point for the MH algorithm.
#' @param npv The number of (effectively returned) plausible values - default is 10.
#' @param approx Whether a normal approximation \code{N(mu,sigma2)} is used to draw the plausible values. Default = TRUE. If FALSE a Metropolitan-Hastings-Algorithm will draw the values.
#' @param thinning A numeric vector of length = 1. If approx = FALSE, a Metropolitan-Hastings-Algorithm draws the plausible values. To avoid autocorrelation, thinning takes every \bold{kth} value as effective plausible value. The default is 6 (every 6th value is taken), which works appropriately in almost all cases here (with default settings).
#' @param burnin How many draws should be discarded at the chains beginning? Default is 10 - and this seems reasonable high (probably 5 will be enough as well), because starting point is the EAP.
#' @param mult Multiplication constant (default = 2). Use this parameter to vary the width of the proposal distribution - which is \code{N(theta_v,mult*SE_eap)} - when a MH-Algorithm is applied. So the constant quantifies the width in terms of multiples of the EAP standard error. 2 works fine with the default thinning. If the supplied value is large, thinning can take lower values without causing autocorrelation.
#' @param ... More arguments
#'
#'
#' @return The function returns a list which main element is \code{pvdraws}. This is a matrix of size number_of_persons x npv - so if 10 plausible values are requested for 100 persons, a 100x10 matrix is returned.
#'
#'@references Mislevy, R. J. (1991). Randomization-based inference about latent variables from complex samples. Psychometrika, 56(2), 177-196.
#'
#'Von Davier, M., Gonzalez, E., & Mislevy, R. (2009). What are plausible values and why are they useful. IERI monograph series, 2, 9-36.
#'
#'Kruschke, J. (2010). Doing Bayesian data analysis: A tutorial introduction with R. Academic Press.
#'
#'
#' @seealso \link{PP_gpcm}, \link{PP_4pl}, \link{JKpp}
#'
#' @rdname PV
#' @example ./R/.example_pv.R
#' @export
#' @author Manuel Reif

PV <- function(estobj,...) UseMethod("PV")


#' @rdname PV
#' @method PV fourpl
#' @export
PV.fourpl <- function(estobj,npv=10,approx=TRUE,thinning=6,burnin=10,mult=2,...)
{
  
  call <- match.call()  
  attr(call, "date") <- date() 
  attr(call,"version") <- packageVersion("PP")
  
  # grep objects 
  respm <- estobj$ipar$respm
  thres <- estobj$ipar$thres
  slopes <- estobj$ipar$slopes
  lowerA <- estobj$ipar$lowerA
  upperA <- estobj$ipar$upperA
  theta_start <- estobj$ipar$theta_start
  mu <- estobj$ipar$mu
  sigma2 <- estobj$ipar$sigma2
  cont <- estobj$ipar$cont
  
  ppresults <- estobj$resPP$resPP

  if(approx)  
  {
    
    pvs <- sapply(1:nrow(ppresults), function(gretel)
    {
      rnorm(npv,ppresults[gretel,1],ppresults[gretel,2])
    })
    
  } else 
  {
    
    # metropolitan-hastings-algorithm  ########
    #################################################
    jederdritte <- (1:(npv*thinning))[1:(npv*thinning) %% thinning == 0]
    
    pvs <- sapply(1:nrow(respm), function(gr)
    {
      
      theta <- ppresults[gr,1] # eap estimate
      # 20 for burnin
      PVvec <- vector(length=npv*thinning+burnin,mode="numeric")
      
      lauf <- 1
      #       zaehl <- 1
      while(lauf <= length(PVvec))
        {
          
          
          Li <- LIK4pl(awv=respm[gr,], thres=thres, slopes=slopes,
                       lowerA=lowerA, upperA=upperA, theta=theta) 
          Post <- Li * dnorm(theta)
          
          # ---- prop >>>
          
          proposed <- rnorm(5,theta,mult*ppresults[gr,2])
          
          for(PROP in proposed)
            {
              
            Li1 <- LIK4pl(awv=respm[gr,], thres=thres, slopes=slopes,
                          lowerA=lowerA, upperA=upperA, theta=PROP) 
            
            Post1 <- Li1 * dnorm(PROP)
            
            # ---- move? >>>
            
            Pmove <- Post1/Post
            
            if(runif(1) <= Pmove)
              {
                theta <- PROP
                PVvec[lauf] <- theta
                lauf <- lauf + 1
                break
              }
  
            }
          
   
        }
        
      PVvec[-(1:burnin)][jederdritte] 
    })
    
    
  }  
  

pvdraws <- t(pvs)  

retpv <- list(pvdraws=pvdraws,call=call)
class(retpv) <- "pv"
return(retpv)

}





#' @rdname PV
#' @method PV gpcm
#' @export
PV.gpcm <- function(estobj,npv=10,approx=TRUE,thinning=6,burnin=10,mult=2,...)
{
  
  call <- match.call()  
  attr(call, "date") <- date() 
  attr(call,"version") <- packageVersion("PP")
  
  # grep objects 
  respm <- estobj$ipar$respm
  thres <- estobj$ipar$thres
  slopes <- estobj$ipar$slopes
  theta_start <- estobj$ipar$theta_start
  mu <- estobj$ipar$mu
  sigma2 <- estobj$ipar$sigma2
  cont <- estobj$ipar$cont
  
  ppresults <- estobj$resPP$resPP
  
  if(approx)  
  {
    
    pvs <- sapply(1:nrow(ppresults), function(gretel)
    {
      rnorm(npv,ppresults[gretel,1],ppresults[gretel,2])
    })
    
  } else 
  {
    
    # metropolitan-hastings-algorithm  ########
    #################################################
    jederdritte <- (1:(npv*thinning))[1:(npv*thinning) %% thinning == 0]
    
    pvs <- sapply(1:nrow(respm), function(gr)
    {
      
      theta <- ppresults[gr,1] # eap estimate
      # 20 for burnin
      PVvec <- vector(length=npv*thinning+burnin,mode="numeric")
      
      lauf <- 1
      #       zaehl <- 1
      while(lauf <= length(PVvec))
      {
        
        Li <- Likgpcm(respm[gr,], thres, slopes,theta) 
        Post <- Li * dnorm(theta)
        
        # ---- prop >>>
        
        proposed <- rnorm(5,theta,mult*ppresults[gr,2])
        
        for(PROP in proposed)
        {
          Li1 <- Likgpcm(respm[gr,], thres, slopes,PROP) 

          
          Post1 <- Li1 * dnorm(PROP)
          
          # ---- move? >>>
          
          Pmove <- Post1/Post
          
          if(runif(1) <= Pmove)
          {
            theta <- PROP
            PVvec[lauf] <- theta
            lauf <- lauf + 1
            break
          }
          
        }
        
        
      }
      
      PVvec[-(1:burnin)][jederdritte] 
    })
    
    
  }  
  
  
  pvdraws <- t(pvs)  
  
  retpv <- list(pvdraws=pvdraws,call=call)
  class(retpv) <- "pv"
  return(retpv)
}





#' @rdname PV
#' @method PV gpcm4pl
#' @export
PV.gpcm4pl <- function(estobj,npv=10,approx=TRUE,thinning=6,burnin=10,mult=2,...)
{
  
  call <- match.call()  
  attr(call, "date") <- date() 
  attr(call,"version") <- packageVersion("PP")
  
  # grep objects 
  respm <- estobj$ipar$respm
  thres <- estobj$ipar$thres
  slopes <- estobj$ipar$slopes
  lowerA <- estobj$ipar$lowerA
  upperA <- estobj$ipar$upperA
  theta_start <- estobj$ipar$theta_start
  mu <- estobj$ipar$mu
  sigma2 <- estobj$ipar$sigma2
  cont <- estobj$ipar$cont
  model2est <- estobj$ipar$model2est
  
  ## where is which model?
  wheregpcm <- model2est == "GPCM"
  where4pl  <- model2est == "4PL"
  
  respm_gpcm <- respm[,wheregpcm, drop=FALSE]
  respm_4pl  <- respm[,where4pl , drop=FALSE]
  
  # 4pl part
  thres4pl  <- thres[1:2,where4pl,drop=FALSE]
  slopes4pl  <- slopes[where4pl]
  lowerA4pl <- lowerA[where4pl]
  upperA4pl <- upperA[where4pl]
  # gpcm part
  thresgpcm <- thres[,wheregpcm,drop=FALSE]
  slopegpcm <- slopes[wheregpcm]
  
  ppresults <- estobj$resPP$resPP
  
  if(approx)  
  {
    
    pvs <- sapply(1:nrow(ppresults), function(gretel)
    {
      rnorm(npv,ppresults[gretel,1],ppresults[gretel,2])
    })
    
  } else 
  {
    
    # metropolitan-hastings-algorithm  ########
    #################################################
    jederdritte <- (1:(npv*thinning))[1:(npv*thinning) %% thinning == 0]
    
    pvs <- sapply(1:nrow(respm), function(gr)
    {
      
      theta <- ppresults[gr,1] # eap estimate
      # 20 for burnin
      PVvec <- vector(length=npv*thinning+burnin,mode="numeric")
      
      lauf <- 1
      #       zaehl <- 1
      while(lauf <= length(PVvec))
      {
        
        Li4pl     <- LIK4pl(awv=respm_4pl[gr,], thres=thres4pl, slopes=slopes4pl,
                             lowerA=lowerA4pl, upperA=upperA4pl, theta=theta) 
        
        Ligpcm    <- Likgpcm(respm_gpcm[gr,],thresgpcm,slopegpcm,theta)
        
        # all together
        Li <- Li4pl * Ligpcm
        
        Post <- Li * dnorm(theta)
        
        # ---- prop >>>
        
        proposed <- rnorm(5,theta,mult*ppresults[gr,2])
        
        for(PROP in proposed)
        {
          
          Li4pl1     <- LIK4pl(awv=respm_4pl[gr,], thres=thres4pl, slopes=slopes4pl,
                              lowerA=lowerA4pl, upperA=upperA4pl, theta=PROP) 

          Ligpcm1    <- Likgpcm(respm_gpcm[gr,],thresgpcm,slopegpcm,PROP)
          
          # all together
          Li1 <- Li4pl1 * Ligpcm1
          
          Post1 <- Li1 * dnorm(PROP)
          
          # ---- move? >>>
          
          Pmove <- Post1/Post
          
          if(runif(1) <= Pmove)
          {
            theta <- PROP
            PVvec[lauf] <- theta
            lauf <- lauf + 1
            break
          }
          
        }
        
        
      }
      
      PVvec[-(1:burnin)][jederdritte] 
    })
    
    
  }  
  

  pvdraws <- t(pvs)  
  retpv <- list(pvdraws=pvdraws,call=call)
  class(retpv) <- "pv"
  return(retpv)
}














########## estimate likelihood of 4pl
LIK4pl <- function(awv, thres, slopes, lowerA, upperA, theta)
{
  
  P <- lowerA + (upperA - lowerA) * exp(slopes*(theta - thres[-1,]))/(1+exp(slopes*(theta - thres[-1,])))
  Q <- 1-P
  
  Li <- P*awv + Q*(1-awv)  
  prod(Li)
  
}





