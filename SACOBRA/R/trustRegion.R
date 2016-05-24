#trustRegion.R
#'
#' Perform trust region refinement
#'
#' If \code{cobra$TrustRegion==TRUE}  (see \code{\link{cobraInit}}), then \code{trustRegion}
#' is called after every iteration in order to refine the best solution so far.
#' This function builds a local model around the best solution and runs a local search in the trust region 
#' to refine the best solution and find a better solution in the neighborhood.
#'
#'
#' @param cobra     an object of class \code{cobra}, which is basically a list  (see \code{\link{cobraInit}})
#' @return the modified \code{cobra} with new/updated elements
#'    \describe{
#'      \item{TRDONE}{ logical, is \code{TRUE} if there are more than d+1 points in the trusted 
#'            region and thus surrogates can be trained. Otherwise \code{FALSE}.}
#'      \item{refinedX}{ if \code{TRDONE==TRUE} the refined solution from the trust-region run,
#'            otherwise \code{NA} }
#'    }
#'    If \code{TRDONE==TRUE} the relevant lists and counters  (\code{A,Fres,df,...}) 
#'    of \code{cobra} will be updated in \code{\link{cobraPhaseII}} as well.
#'    
#' @author Samineh Bagheri (\email{samineh.bagheri@@fh-koeln.de})

trustRegion<-function(cobra){
  #cobra$radi is calculated according to the progress in the last iterations by adaptRadi function
  delta<-cobra$radi[length(cobra$radi)]
  TRlower<-pmax(cobra$xbest-delta,cobra$lower)
  TRupper<-pmin(cobra$xbest+delta,cobra$upper)
  
 
   #indices of the points in the trust region hypercube
  ind<-which(sapply(1:nrow(cobra$A),indCheck<-function(i){
    (TRupper-cobra$A[i,] >= 0 ) && (cobra$A[i,]-TRlower >= 0)
  }))
  
   dimension<-ncol(cobra$A)
   if(length(ind) < dimension+1 ){ # for building an RBF we need at least d+1 points, if we have less points then building 
      refinedX<-NA            # a local model in the defined trust region is not possible
      print("Trust Region cannot be performed ")
      cobra$TRDONE<-FALSE
  }else{ 
    cobra$TRDONE<-TRUE
    #build the local model over the trust region
  }
    
  if(cobra$TRDONE){
    ############################Training the local model###############################################
    cat(paste(">> Training Trust Region" ,cobra$RBFmodel,"surrogates","...\n" ))
    TRA<-cobra$A[ind,]
    TRGres<-cobra$Gres[ind,]
    TRFres<-cobra$Fres[ind]
    # TRFres<-sapply(1:length(TRFres) , function(i){rescale(TRFres[i],to=c(-1,1),from=c(min(TRFres),max(TRFres)))
    # })
    #TRFres<-TRFres/max(TRFres)
    sw=switch(cobra$RBFmodel,
              "cubic" =     {constraintSurrogatesTR <- trainCubicRBF(TRA,TRGres,squares=cobra$squaresC)
                             fitnessSurrogateTR <- trainCubicRBF(TRA,TRFres,squares=cobra$squaresF)},
              "Gaussian"=   {constraintSurrogatesTR <- trainGaussRBF(TRA,TRGres,cobra$RBFwidth,squares=cobra$squaresC);
                             fitnessSurrogateTR <- trainGaussRBF(TRA,TRFres,cobra$RBFwidth,squares=cobra$squaresF)}
    )
    
    #trainTRsurrogates
    ############################Training the local model finished###############################################
    
    gama<-0
    ro<-gama*cobra$l  
    
    
    
    subProbTR2 <- function(x){
      y<-predict.RBFinter(fitnessSurrogateTR,matrix(x,ncol=dimension))
      return(y)
    }
    
    gCOBRATR <- function(x) {
      h <- c()
      distance <- distLine(x,cobra$A[ind,])
      subC<-pmax((ro-distance),0)
      h[1] <- sum(subC)*cobra$drFactor
      #cat(">>> Entering interpRBF\n")
      constraintPredictionTR <- interpRBF(x,constraintSurrogatesTR)+cobra$epsilonInit^2
      h <- (-1.0)*c(h[1], constraintPredictionTR) # TODO -1* ... is required for COBYLA constraints, maybe also for other optimizers? 
      #cat("<<< leaving interpRBF\n")
      return(h)
    }
    
    #generate a random starting point in the trust region
    xStart<-runif(ncol(cobra$A),min=TRlower,max=TRupper)
    # xStart<-cobra$xbest
    switch(cobra$seqOptimizer,
           COBYLA={ subMin<-nloptr::cobyla(xStart,fn=subProbTR2,lower=TRlower,upper=TRupper,hin=gCOBRATR,control=list(maxeval=cobra$seqFeval,xtol_rel=cobra$seqTol)); subMin$feval=subMin$iter },
           ISRES ={ subMin<-isres2(xStart,fn=subProbTR2,lower=TRlower,upper=TRupper,hin=gCOBRATR, maxeval=cobra$seqFeval); subMin$feval=subMin$iter}
           #RANDOMSEARCH = { subMin <- randomSearch(xStart, fn=subProb2, lower=TRlower,upper=TRupper,control=list(maxfeval=cobra$seqFeval, sd=0.05))}
    )
    #browser()
    refinedX<-subMin$par
  }
  cobra$refinedX<-refinedX

  return(cobra)
}

