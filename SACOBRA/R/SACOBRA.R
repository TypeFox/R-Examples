#
#
#Samineh Bagheri 27.04.2015
#
#
#SACOBRA.R
#
#
#

######################################################################################
######################################################################################
# Package Description for Roxygene:
#' Self-adjusting Constrained Optimization with RBF 
#'
#' \tabular{ll}{
#' Package: \tab SACOBRA\cr
#' Type: \tab Package\cr
#' Version: \tab 0.7\cr
#' Date: \tab 30.08.2015\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' SACOBRA is a package for numeric constrained optimization of expensive black-box functions under severely 
#' limited budgets. It is an extension of the COBRA algorithm by Regis (R. Regis: "Constrained 
#' optimization by radial basis function interpolation for high-dimensional expensive black-box 
#' problems with infeasible initial points", Engineering Optimization, Taylor & Francis, 46, p. 218-243, 2013)
#' 
#' These extensions include: \cr
#' 1) A repair algorithm for infeasible solutions, \cr
#' 2) several internal optimizers and several initial design generation methods,  \cr
#' 3) self-adjusting random restart algorithm,  \cr
#' 4) self-adjusting logarithmic transform for objective functions with large output ranges,  \cr
#' 5) range normalization of constraint functions, \cr
#' 6) self-adjusting DRC selection. \cr
#' 
#' SACOBRA performs optimization with a minimum of true function evaluations. It has proven
#' to work well on probelms with high dimensions (e.g. d=124) and many constraints (e.g. 60).
#' It is usable for all kind of numeric optimization, but not for combinatorial optimization.
#' 
#' For more details see:\cr
#' Koch, P.; Bagheri, S.; Konen, W. et al.: "A New Repair Method For Constrained Optimization". 
#' In: Proceedings of the 17th Genetic and Evolutionary Computation Conference, 2015, 
#' \url{http://www.gm.fh-koeln.de/~konen/Publikationen/Koch2015a-GECCO.pdf} \cr
#' and \cr
#' Koch, P.; Bagheri, S. et al.: "Constrained Optimization with a Limited Number of Function Evaluations"
#' In: W. Hoffmann, F. & Huellermeier, E. (Eds.),  Proceedings 24. Workshop Computational Intelligence, 
#' Universitaetsverlag Karlsruhe, 2014, 119-134, \url{http://www.gm.fh-koeln.de/~konen/Publikationen/Koch2014a-GMA-CI.pdf}.
#' 
#' The main entry point functions are \code{\link{cobraInit}} and \code{\link{startCobra}}. 
#' See \code{\link{startCobra}} for an example and \code{\link{cobraInit}} for an overview of adjustable SACOBRA-parameters.
#' Another example is in \code{\link{multiCOBRA}}.
#'                                                
#' @name SACOBRA-package
#' @aliases SACOBRA
#' @docType package
#' @title Self-adjusting Constrained Optimization with RBF
#' @author Samineh Bagheri (\email{Samineh.Bagheri@@fh-koeln.de}), Wolfgang Konen (\email{Wolfgang.Konen@@fh-koeln.de}), Patrick Koch
#' @references \url{http://lwibs01.gm.fh-koeln.de/blogs/ciop/research/monrep/}
#' @keywords package constraints optimization RBF
#' @import testit
#' @import grDevices
#' @import graphics
#' @import methods
#' @import stats
#' @import utils

#End of Package Description
NA #NULL, ends description without hiding first function
######################################################################################
######################################################################################


# long and short DRC
DRCL   <-c(0.3,0.05, 0.001, 0.0005,0.0)
DRCS   <-c(0.001,0.0)

#   rescale borders ### deprecated, use cobra$newlower, cobra$newupper
#LOWER  <- -1.0    # WK: Bug Fix, was 0.0
#UPPER  <- 1.0


#
rescaleWrapper<-function(fn,lower,upper,dimension,newlower,newupper){
  oldfn<-fn                     
  newfn<-function(x){
    x<-sapply(1:length(x) , function(i){scales::rescale(x[i],to=c(lower[i],upper[i]),from=c(newlower,newupper))
    })
    y<-oldfn(x)
    return(y)
  }
  return(newfn)
}

#' rescaling
#' 
#' Scale vector x in original space forward to rescaled space (usually \eqn{[-1,1]^d})
#' 
#' @param x       a vector in the original input space 
#' @param cobra   list from \code{\link{cobraInit}}, we need here
#' \describe{
#'    \item{originalL}{   a vector with lower bounds in original input space}
#'    \item{originalU}{   a vector with upper bounds in original input space}
#'    \item{newlower}{    a number, the rescaled lower bound for all dimensions}
#'    \item{newupper}{    a number, the rescaled upper bound for all dimensions}
#' }
#' 
#' @return \code{z},      scaled version of vector x
#' @seealso   \code{\link{inverseRescale}}
#' @export
#' 
forwardRescale <- function(x,cobra) {
  lb<-rep(cobra$newlower,length(x))
  up<-rep(cobra$newupper,length(x))
  origL <- cobra$originalL
  origU <- cobra$originalU
  z<-sapply(1:length(x) , function(i){scales::rescale(x[i],from=c(origL[i],origU[i]),to=c(lb[i],up[i]))})
  return(z)
}

#' inverse rescaling
#' 
#' Scale vector x in rescaled space back to original space
#' 
#' @param x       a vector in the rescaled input space (usually \eqn{[-1,1]^d})
#' @param cobra   list from \code{\link{cobraInit}}, we need here
#' \describe{
#'    \item{originalL}{   a vector with lower bounds in original input space}
#'    \item{originalU}{   a vector with upper bounds in original input space}
#'    \item{newlower}{    a number, the rescaled lower bound for all dimensions}
#'    \item{newupper}{    a number, the rescaled upper bound for all dimensions}
#' }
#' 
#' @return \code{z},      inverse rescaling of vector x
#' @seealso   \code{\link{forwardRescale}}
#' @export
#' 
inverseRescale <- function(x,cobra) {
  lb<-rep(cobra$newlower,length(x))
  up<-rep(cobra$newupper,length(x))
  z<-sapply(1:length(x) , function(i){scales::rescale(x[i],to=c(cobra$originalL[i],cobra$originalU[i]),from=c(lb[i],up[i]))})
  return(z)
}

#' Return best feasible solution in original space
#' 
#' @param cobra an object of class COBRA (see  \code{\link{cobraInit}})
#' 
#' @return the best feasible solution in original space
#' @seealso   \code{\link{getFbest}}
#' @export
#' 
getXbest <- function(cobra) {
  xbest <- cobra$xbest
  if (cobra$rescale) xbest <- inverseRescale(xbest,cobra)
  return(xbest)
}

#' Return best objective function value
#'
#' Return the original objective function value at the best feasible solution 
#' 
#' Note: We cannot take the best function value via \code{cobra$fn}, because this 
#' may be modified by \code{plog()} or others )
#' 
#' @param cobra an object of class COBRA (see  \code{\link{cobraInit}})
#' 
#' @return the original objective function value at the best feasible solution 
#' @seealso   \code{\link{getXbest}}
#' @export
#' 
getFbest <- function(cobra) {
  return(cobra$originalfn(getXbest(cobra))[1])
}

#Random start Algorithm
#
RandomStart<-function(cobra){
  anewrand<-runif(1,min=0,max=1)
  
 # randomnessTemp<-(9*0.3/20)*tanh(-(nrow(cobra$A)-(cobra$initDesPoints+15)))+11*0.3/20
 diff<-cobra$sac$RSmax-cobra$sac$RSmin
 switch(cobra$sac$RStype,
        SIGMOID =  randomnessTemp<-(diff/2)*tanh(-(nrow(cobra$A)-(cobra$initDesPoints+15)))+(diff/2)+cobra$sac$RSmin,
        CONSTANT= randomnessTemp<-diff/2
)

  #cat(paste(randomnessTemp,"=randomness"))

  if( (anewrand< randomnessTemp)  || (cobra$progressCount >= cobra$sac$Cs)){
    verboseprint(cobra$verbose, important=FALSE,"Starting the internal optimizer with a random point in the space")
    xStart<-runif(n=length(cobra$xbest),min=cobra$lower,max=cobra$upper)
    #browser()
   cobra$progressCount<-0
  } else{
    xStart<-cobra$xbest
  }
  cobra$xStart<-xStart
  return(cobra)
}


#Adjust DRC
#
adDRC<-function(maxF,minF){
  FRange<-(maxF-minF)
  if(FRange>1e+03) {
    cat(sprintf("FR=%g is large, XI is set to Short DRC \n",FRange))
    DRC<-DRCS
  }else{
    DRC<-DRCL
    cat(sprintf("FR=%g is not large, XI is set to Long DRC\n",FRange))
    
  }
  return(DRC)
}

#Adjust Constraint functions
#
adCon<-function(cobra){
  detLen <-function(x){
    maxL<-max(x)
    minL<-min(x)
    return(maxL-minL)
  }
  fnold<-cobra$fn
  
  
  GRL<-apply(cobra$Gres,2,detLen)
  if (min(GRL)==0) {  # pathological case where at least one constraint is constant:
    GR <- -Inf        # inhibit constraint normalization
  } else {
    GR<-max(GRL)/min(GRL)
  }
  if(GR > cobra$sac$TGR){
    cat(sprintf("Normalizing Constraint Functions \n"))
    GRfact<-c(1,GRL*(1/mean(GRL)))
    fn<-function(x){
      return(fnold(x)/GRfact)
    }
    cobra$fn<-fn

    Gres<-NULL
    for(i in 1:nrow(cobra$Gres)){
      Gres<-rbind(Gres,cobra$Gres[i,]/GRfact[-1])
    }
    cobra$Gres<-Gres
    
  }
  
  return(cobra)
}

#Adjust Fitness Function
#
# Note that the fitness function cobra$fn is not changed by this function. Instead, all results
# found in cobra$Fres are transformed with transfFunc and the transformed results are saved 
# on cobra$SurrogateInput. This is then used to train a new surrogate model for the fitness 
# function and use this model in the sequential optimizer.
#
adFit<-function(cobra,ind){
  maxF=max(cobra$Fres)
  minF=min(cobra$Fres)
  FRange<-(maxF-minF)
  
  transfFunc<-function(x,pShift=0){
    y<-plog(x,pShift=pShift)
    return(y)
  }
  
  #If cobra$online PLOG is true then the desciosn to do the plog transfomation or not 
  #is being made in every iteration according to the p-effect otherwise the decision is made once accoridng to the FRange value 
  if(cobra$sac$onlinePLOG){
    pShift<-0
    #decision making according to p-effect
    if(cobra$pEffect > 1){
      Fres<-sapply(cobra$Fres[ind],transfFunc,pShift=pShift)
      cobra$PLOG<-c(cobra$PLOG,TRUE)
    }else{
      Fres<-cobra$Fres[ind]
      pShift<-NA
      cobra$PLOG<-c(cobra$PLOG,FALSE)
    }
  }else{
    if(FRange>cobra$sac$TFRange) {   
      if(cobra$sac$adaptivePLOG){
        # browser()
        verbosecat(verbose=cobra$verbose,important=cobra$important,sprintf("Very Large FR=%g, applying adaptive plog()",FRange))
        #pShift<-mean(c(cobra$fbest,min(cobra$Fres)[1]))
        pShift<-mean(c(cobra$fbest,0))   # /SB/just for testing purposes
        #pShift<-cobra$fbest               # /WK/ the most simple alternative (don't know quality yet)
        
        #transfFunc<-function(x,pShift){
        #      x<-x-mean(c(pShift,0))
        #      #The following print is only for debugging purposes
        #      #print(paste("pShitf:",pShift))
        #      y<-plog(x,pShift=0)
        #     return(y)
        #    }
      }else{
        verbosecat(cobra$verbose,sprintf("Very Large FR=%g, applying plog()\n",FRange))
        pShift<- 0
      }
      Fres<-sapply(cobra$Fres[ind],transfFunc,pShift=pShift)
      cobra$PLOG<-TRUE
      
    }else{
      Fres<-cobra$Fres[ind]
      pShift<-NA
      cobra$PLOG<-FALSE
      
    } 
  }
  
  cobra$pShift<-c(cobra$pShift,pShift)
  cobra$SurrogateInput<-Fres
  return(cobra)
}

#' Monotonic transform
#' 
#' The function is introduced in [Regis 2014] and extended here by a parameter \eqn{p_{shift}}.\cr 
#' Let \eqn{y' = (y-p_{shift})}: 
#'  \deqn{ plog(y) =  \ln(1+ y'), \quad\mbox{if}\quad y' \geq 0 } 
#'  \deqn{ plog(y) = -\ln(1- y'), \quad\mbox{if}\quad y'  <   0 } 
#'  
#' @param y       function argument
#' @param pShift  shift
#'  
#' @return \eqn{plog(y)}
#' @seealso \code{\link{plogReverse}}
#' @export
#'  
plog<-function(y,pShift=0.0){
  if(y-pShift>=0){
    ret<- log(1+y-pShift)
  }else{
    ret<- -log(1-(y-pShift))
  }
  #return(100*ret)
  return(ret)
}

#' Inverse of \code{\link{plog}}
#' 
#' @param y       function argument
#' @param pShift  shift
#' 
#' @return \eqn{plog^{-1}(y)}
#' @seealso \code{\link{plog}}
#' @export
plogReverse<-function(y,pShift=0){
  #y <- y/100
  if(y > 0){
    ret<-exp(y)-1+pShift
  }else{
    ret<-pShift+1-(1/exp(y)) 
  } 
  return(ret)
}

# --- alternative plog with sqrt ----
#
# plog<-function(y,pShift=0.0){
#   if(y-pShift>=0){
#     ret<- sqrt(1/4+y-pShift)-sqrt(1/4)
#   }else{
#     ret<- -sqrt(1/4-(y-pShift))+sqrt(1/4)
#   }
#   #return(100*ret)
#   return(ret)
# }
# 
# plogReverse<-function(y,pShift=0){
#   #y <- y/100
#   if(y > 0){
#     ret <- (y+1/2)^2-1/4+pShift
#   }else{
#     ret <- -(y-1/2)^2+1/4+pShift 
#   } 
#   return(ret)
# }
