#' Run a jackknife
#' 
#' This function uses a jackknife approach to compute person parameters. The jackknife ability measure is based on primarily estimated models (\code{PP_4pl()}, \code{PP_gpcm()} or \code{PPall()}) - so the function is applied on the estimation objects, and jackknifed ability measures are returned.
#' 
#' Please use the Jackknife \bold{Standard-Error} output with \bold{caution}! It is implemented as suggested in Wainer and Wright (1980), but the results seem a bit strange, because the JK-SE is supposed to overestimate the SE compared to the MLE-SE. Actually, in all examples an underestimation of the SE was observed compared to the MLE/WLE-SE! 
#' 
#' 
#' \bold{AMT-robustified jackknife:} When choosing \code{cmeth = AMT}, the jackknife ability subsample estimates and the original supplied ability estimate are combined to a single jackknife-ability value by the Sine M-estimator. The AMT (or Sine M-estimator) is one of the winners in the Princeton Robustness Study of 1972. To get a better idea how the estimation process works, take a closer look to the paper which is mentioned below (Wainer & Wright, 1980).
#' 
#' @param estobj An object which originates from using \code{PP_gpcm()}, \code{PP_4pl()} or \code{PPall()}.
#' 
#' @param ...  More input.
#' @export
#' @rdname Jkpp
#' @references
#' Wainer, H., & Wright, B. D. (1980). Robust estimation of ability in the Rasch model. Psychometrika, 45(3), 373-391.
#'
#'
#' @seealso \link{PP_gpcm}, \link{PP_4pl}, \link{PPall}
#' 
#'@example ./R/.examples_JK.R
#'@author Manuel Reif
JKpp <- function(estobj,...) UseMethod("JKpp")


# ---------------------------------------------------------------------


#' @rdname Jkpp
#' @param cmeth Choose the centering method, to summarize the n jackknife results to one single ability estimate. There are three valid entries: "mean", "median" and "AMT" (see Details for further description).
#' @param maxsteps The maximum number of steps the NR Algorithm will take.
#' @param exac How accurate are the estimates supposed to be? Default is 0.001.
#' @param fullmat Default = FALSE. If TRUE, the function returns the whole jackknife matrix, which is the basis for the jackknife estimator.
#' @param ctrl More controls
#' 
#' @method JKpp fourpl
#' @export
JKpp.fourpl <- function(estobj, cmeth="mean", maxsteps=500,
                        exac=0.001, fullmat=FALSE,ctrl=list(), ...)
{
  
call <- match.call()  
attr(call, "date") <- date() 
attr(call,"version") <- packageVersion("PP")

# pick out objects  
respm <- estobj$ipar$respm
thres <- estobj$ipar$thres
slopes <- estobj$ipar$slopes
lowerA <- estobj$ipar$lowerA
upperA <- estobj$ipar$upperA
#theta_start <- estobj$ipar$theta_start
theta_start <- rep(0,nrow(respm))
mu <- estobj$ipar$mu
sigma2 <- estobj$ipar$sigma2
cont <- estobj$ipar$cont
H <- estobj$ipar$H

# in case the killdupli option = TRUE 
## if not, create a vector which indicates merely each row
if(!is.null(estobj$ipar$dupvec$posvec))
  {
    
  POS <- estobj$ipar$dupvec$posvec
  } else 
    {
      POS <- 1:nrow(respm)  
    }
    
type <- estobj$type 
loa <- 1:ncol(respm)
# save the responses into this matrix
jk_mat <- matrix(0,nrow=nrow(respm),ncol=length(loa))

#browser()
  
# run the 4pl jackknife 
if(type %in% c("mle","wle","map","robust"))
  {
    
  for(jkrun in loa)
    {
      
    
    jk_mat[,jkrun] <- NR_4PL(respm[,-jkrun,drop=FALSE],DELTA = thres[,-jkrun,drop=FALSE],ALPHA = slopes[-jkrun], LOWA = lowerA[-jkrun],UPPA = upperA[-jkrun], THETA = theta_start, wm=type,maxsteps,exac,mu,sigma2,H=H)$resPP[,1] 
        
    }
  
  } else if(type == "eap")
      {
    
    
      for(jkrun in loa)
        {
            jk_mat[,jkrun] <- eap_4pl(respm[,-jkrun,drop=FALSE], thres[,-jkrun,drop=FALSE],
                                      slopes[-jkrun], lowerA=lowerA[-jkrun],
                                      upperA=upperA[-jkrun],mu = mu, sigma2 = sigma2)[,1]
        }
      }


notna <- !is.na(estobj$resPP$resPP[,2])

RES <- estobj$resPP$resPP[notna,1]

psvalues  <- RES *ncol(respm) - jk_mat[POS,] * (ncol(respm) - 1)

output_jk <- cco(psvalues,cmeth) 

### jk standard errors - see page 388
L <- ncol(respm)
Lm1 <- L - 1

jk_se <- sqrt(rowSums((jk_mat[POS,] - output_jk)^2,na.rm=TRUE)/(L*Lm1))


resjk_raw <- cbind(output_jk,jk_se)

if(cont$killdupli)
  {
    resjk_raw <- resjk_raw[POS,]
  }
  

resjk <- estobj$resPP$resPP
resjk[notna,] <- resjk_raw
  
if(fullmat)
  {
  
  endres <- list(resjk=resjk,call=call,type=type,jk_mat=jk_mat)  
    
  } else 
    {
      
    endres <- list(resjk=resjk,call=call,type=type)
      
    }

class(endres) <- c("jk")
return(endres)

}







# ---------------------------------------------------------------------

#' @rdname Jkpp
#' @method JKpp gpcm
#' @export
JKpp.gpcm <- function(estobj, cmeth="mean", maxsteps=500, 
                      exac=0.001, fullmat=FALSE,  ctrl=list(), ...)
{
  
  call <- match.call()  
  attr(call, "date") <- date() 
  attr(call,"version") <- packageVersion("PP")
  
  
  # pick out objects  
  respm <- estobj$ipar$respm
  thres <- estobj$ipar$thres
  slopes <- estobj$ipar$slopes
#   lowerA <- estobj$ipar$lowerA
#   upperA <- estobj$ipar$upperA
  #theta_start <- estobj$ipar$theta_start
  theta_start <- rep(0,nrow(respm))
  mu <- estobj$ipar$mu
  sigma2 <- estobj$ipar$sigma2
  cont <- estobj$ipar$cont
  H <- estobj$ipar$H  


if(!is.null(estobj$ipar$dupvec$posvec))
  {
    
    POS <- estobj$ipar$dupvec$posvec
  } else 
    {
      POS <- 1:nrow(respm)  
    }

# skip -inf   

  type <- estobj$type 
  loa <- 1:ncol(respm)
  # save the responses into this matrix
  jk_mat <- matrix(0,nrow=nrow(respm),ncol=length(loa))
  
  


if(type %in% c("mle","wle","map","robust"))
{
  
  for(jkrun in loa)
  {
    
    
    jk_mat[,jkrun] <- NR_GPCM(respm[,-jkrun,drop=FALSE], thres[,-jkrun,drop=FALSE],
                              slopes[-jkrun], theta_start, type,
                              maxsteps, exac, mu, sigma2,H=H)$resPP[,1] 
    
  }
  
} else if(type == "eap")
    {
      
      for(jkrun in loa)
        {
        jk_mat[,jkrun] <- eap_gpcm(respm[,-jkrun,drop=FALSE], thres[,-jkrun,drop=FALSE], slopes[-jkrun],
                                   mu = mu, sigma2 = sigma2)[,1]
        }
    }





  
  notna <- !is.na(estobj$resPP$resPP[,2])
  
  RES <- estobj$resPP$resPP[notna,1]
  
  psvalues  <- RES *ncol(respm) - jk_mat[POS,] * (ncol(respm) - 1)
  
  output_jk <- cco(psvalues,cmeth) 
  
  ### jk standard errors - see page 388
  L <- ncol(respm)
  Lm1 <- L - 1
  
  jk_se <- sqrt(rowSums((jk_mat[POS,] - output_jk)^2,na.rm=TRUE)/(L*Lm1))
  
  resjk_raw <- cbind(output_jk,jk_se)
  
  if(cont$killdupli)
    {
      resjk_raw <- resjk_raw[POS,]
    }

  
  resjk <- estobj$resPP$resPP
  resjk[notna,] <- cbind(output_jk,jk_se)
  
  if(fullmat)
    {
      
      endres <- list(resjk=resjk,call=call,type=type,jk_mat=jk_mat)  
      
    } else 
      {
        
        endres <- list(resjk=resjk,type=type,call=call)
        
      }

  class(endres) <- c("jk")
  return(endres)

}





# ---------------------------------------------------------------------

#' @rdname Jkpp
#' @method JKpp gpcm4pl
#' @export 
JKpp.gpcm4pl <- function(estobj, cmeth="mean", maxsteps=500,
                         exac=0.001, fullmat=FALSE, ctrl=list(), ...)
{
  
  call <- match.call()  
  attr(call, "date") <- date() 
  attr(call,"version") <- packageVersion("PP")
  
  # pick out objects  
  respm  <- estobj$ipar$respm
  thres  <- estobj$ipar$thres
  slopes <- estobj$ipar$slopes
  lowerA <- estobj$ipar$lowerA
  upperA <- estobj$ipar$upperA
  #theta_start <- estobj$ipar$theta_start
  theta_start <- rep(0,nrow(respm))
  mu <- estobj$ipar$mu
  sigma2 <- estobj$ipar$sigma2
  cont <- estobj$ipar$cont
  H <- estobj$ipar$H
  model2est <- estobj$ipar$model2est

  if(!is.null(estobj$ipar$dupvec$posvec))
    {
      
      POS <- estobj$ipar$dupvec$posvec
    } else 
      {
        POS <- 1:nrow(respm)  
      }
  
  
  type <- estobj$type 
  loa <- 1:ncol(respm)
  # save the responses into this matrix
  jk_mat <- matrix(0,nrow=nrow(respm),ncol=length(loa))
  
  
  
  if(type %in% c("mle","wle","map","robust"))
  {
    
    for(jkrun in loa)
    {
      
      
      jk_mat[,jkrun] <- NR_mixed(awm=respm[,-jkrun,drop=FALSE],DELTA = thres[,-jkrun,drop=FALSE],ALPHA = slopes[-jkrun],LOWA = lowerA[-jkrun],UPPA = upperA[-jkrun],THETA = theta_start, model=model2est,wm=type,maxsteps=maxsteps,
                                 exac=exac,mu=mu,sigma2=sigma2,H=H)$resPP[,1] 
      
    }
    
  } else if(type == "eap")
    {
      
      for(jkrun in loa)
        {
        jk_mat[,jkrun] <- eap_mixed(respm[,-jkrun,drop=FALSE], thres[,-jkrun,drop=FALSE],
                                    slopes[-jkrun], lowerA=lowerA[-jkrun],
                                    upperA=upperA[,-jkrun], mu = mu, sigma2 = sigma2,
                                    model2est=model2est)[,1]
        }
    }
    
  

  

  notna <- !is.na(estobj$resPP$resPP[,2])
  
  RES <- estobj$resPP$resPP[notna,1]
  
  psvalues  <- RES *ncol(respm) - jk_mat[POS,] * (ncol(respm) - 1)
  
  output_jk <- cco(psvalues,cmeth) 
  
  ### jk standard errors - see page 388
  L <- ncol(respm)
  Lm1 <- L - 1
  
  jk_se <- sqrt(rowSums((jk_mat[POS,] - output_jk)^2,na.rm=TRUE)/(L*Lm1))
  
  resjk_raw <- cbind(output_jk,jk_se)
  
  if(cont$killdupli)
    {
      resjk_raw <- resjk_raw[POS,]
    }
    
  
  resjk <- estobj$resPP$resPP
  resjk[notna,] <- cbind(output_jk,jk_se)
  
  
  if(fullmat)
    {
      
      endres <- list(resjk=resjk,call=call,type=type,jk_mat=jk_mat)  
      
    } else 
      {
        
        endres <- list(resjk=resjk,type=type,call=call)
        
      }
  class(endres) <- c("jk")
  return(endres)
  
}










####### compute AMT
AMTnew <- function(TT,xj)
{
  suppressWarnings(zw <- sin((xj + TT)/2.1))
  nullk <- abs(xj) < 2.1*pi
  zw[!nullk] <- 0
  
  sum(zw)
}


###### compute composed PP

cco <- function(psvalues,cmeth)
{

  if(cmeth=="mean")
  {
    
    jkest <- rowMeans(psvalues,na.rm=TRUE)
    
  } else if(cmeth=="median")
  {
    jkest <- apply(psvalues,1,function(cen) median(cen,na.rm=TRUE))
    
  } else if(cmeth=="AMT")
  {
    
    jkest <- apply(psvalues,1,function(psv)
    {
      negpv <- psv[psv <= 0 & !is.na(psv) & !is.nan(psv)]
      pospv <- psv[psv > 0 & !is.na(psv) & !is.nan(psv)]
      
      if(all(is.nan(psv)))
      {
        yes <- NaN
      } else 
      {
        # wir haben hier bei MLE noch ein problem wenn NaN in psvalues enthalten sind. dasselbe problem gilt auch fuer wle, nur wird das hier wesentlich seltener vorkommen.
        an1 <- optim(0,AMTnew,xj=pospv,method="BFGS",control = list(fnscale=-1)) # minimize
        an2 <- optim(0,AMTnew,xj=negpv,method="BFGS",control = list(fnscale=1)) # maximize
        
        yes <- sum(an1$par*length(pospv) + an2$par*length(negpv))/length(psv)
        yes
        
      }
      
    })
    
  }    
  
return(jkest)  
}












