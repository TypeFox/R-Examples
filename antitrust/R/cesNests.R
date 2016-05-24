

setClass(
         Class   = "CESNests",
         contains="CES",

         representation=representation(
         nests="factor",
         parmsStart="numeric",
         constraint="logical"
         ),

         prototype=prototype(
         parmsStart      =  numeric(),
         constraint=TRUE
         ),

         validity=function(object){




             nprods    <- length(object@prices)
             nNestParm <- nlevels(object@nests) #calculate the number of nesting parameters
             nMargins  <- length(object@margins[!is.na(object@margins)])
             maxNests  <- nMargins - 1

             ## Identify Singleton Nests
             nestCnt   <- tapply(object@prices,object@nests,length)
             nestCnt   <- nestCnt[object@nests]
             isSingleton <- nestCnt==1

             nNestParm <- nNestParm - sum(isSingleton) #singleton nests are not identified

             if(identical(nNestParm,1)) stop("'ces.nests' cannot be used for non-nested problems or problems with only singleton nests. Use 'ces' instead")

             if(nprods != length(object@nests)){
                 stop("'nests' length must equal the number of products")
             }

              if(object@constraint && length(object@parmsStart)!=2){
                 stop("when 'constraint' is TRUE, 'parmsStart' must be a vector of length 2")
                 }
             else if(!object@constraint && nNestParm + 1 != length(object@parmsStart)){
                 stop("when 'constraint' is FALSE, 'parmsStart' must be a vector of length ",nNestParm + 1)
             }


             if(!object@constraint &&
                any(tapply(object@margins[!isSingleton],object@nests[!isSingleton],
                           function(x){if(all(is.na(x))){return(TRUE)} else{return(FALSE)}}
                           )
                    ,na.rm=TRUE)
                ){
                 stop("when 'constraint' is FALSE, at least one product margin must be supplied for each non-singleton nest")
             }


             if(nNestParm > nMargins){
                 stop(paste(
                            "Impossible to calibrate nest parameters with the number of margins supplied.\n",
                            "The maximum number of nests supported by the supplied margin information is"
                            ,maxNests,"."))
             }
         }

         )


setMethod(
          f= "calcSlopes",
          signature= "CESNests",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex
              shareInside  <-  object@shareInside
              nests        <- object@nests
              parmsStart   <- object@parmsStart
              constraint   <- object@constraint

              nestCnt      <- tapply(prices,nests,length)


              isSingletonNest <- nestCnt==1

              if(any(isSingletonNest)){
                  warning("Some nests contain only one product; their nesting parameters are not identified.
Normalizing these parameters to 1.")

              }



              if(!constraint){
                  parmsStart   <- parmsStart[c(TRUE,!isSingletonNest)] #always retain first element; this is
                                                                          # the initial value for price coefficient
              }

               ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              ## identify which products have enough margin
              ## information to impute Bertrand margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))

              sharesNests <- tapply(shares,nests,sum)[nests]

              sharesNests <- shares / sharesNests


              ## back out the parameter on the numeraire, when appropriate
              if(shareInside<1) {alpha <- 1/shareInside -1}
              else{ alpha <- NULL}






              ## Estimate parameters by
              ## Minimizing the distance between observed and predicted margins
              minD <- function(theta){

                  gamma  <- theta[1]

                  sigma <- as.numeric(!isSingletonNest) # normalize singleton nest parms to 0
                  sigma[!isSingletonNest] <- theta[-1]

                  elast <- diag(sigma - gamma)

                  elast <- elast[nests,nests]
                  elast <- elast * matrix(sharesNests,ncol=nprods,nrow=nprods)
                  elast <- elast + (gamma-1) * matrix(shares,ncol=nprods,nrow=nprods)

                  diag(elast) <- diag(elast) - sigma[nests]

                  #marginsCand <- -1 * as.vector(ginv(elast * ownerPre) %*% (shares * diag(ownerPre))) / shares
                  #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)


                  elast      <- elast[isMargin,isMargin]
                  shares     <- shares[isMargin]
                  ownerPre   <- ownerPre[isMargin,isMargin]
                  margins    <- margins[isMargin]

                  FOC <- (shares * diag(ownerPre)) + (elast * ownerPre) %*% (shares * margins)
                  measure<-sum(FOC^2,na.rm=TRUE)

                  return(measure)
              }

              ## Constrain optimizer to look for solutions where sigma_i > gamma > 1 for all i
              constrA <- diag(length(parmsStart))
              constrA[-1,1] <- -1

              constrB <- rep(0,length(parmsStart))
              constrB[1] <- 1

              minTheta <- constrOptim(parmsStart,minD,grad=NULL,ui=constrA,ci=constrB)


              if(minTheta$convergence != 0){
                  warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
              }





              minGamma <- minTheta$par[1]
              names(minGamma) <- "Gamma"

              minSigma <-  as.numeric(!isSingletonNest)
              minSigma[!isSingletonNest] <- minTheta$par[-1]


              minSigmaOut        <- minSigma
              minSigma           <- minSigma[nests]
              names(minSigmaOut)    <- levels(nests)


               if(is.na(idx)){
                  idxShare      <- 1 - sum(shares)
                  idxShareNests <- 1
                  idxPrice      <- object@priceOutside
                  idxSigma      <- 0
              }

              else{

                  idxShare      <- shares[idx]
                  idxShareNests <- sharesNests[idx]
                  idxPrice      <- prices[idx]
                  idxSigma      <- minSigma[idx]
               }


              meanval <-
                  log(shares) - log(idxShare) + (minGamma - 1) *
                  (log(prices) - log(idxPrice)) -
                      (minSigma-minGamma)/(minSigma-1)*log(sharesNests) +
                         (idxSigma - minGamma)/(idxSigma-1)*log(idxShareNests)

              meanval <- exp( (minSigma-1)/(minGamma-1) * meanval )


              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=alpha,gamma=minGamma,sigma=minSigmaOut,meanval=meanval)

              return(object)
          }
          )


setMethod(
 f= "calcShares",
 signature= "CESNests",
 definition=function(object,preMerger=TRUE,revenue=FALSE){



     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}


     nests     <- object@nests
     gamma    <- object@slopes$gamma
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval

     outVal <- ifelse(sum(object@shares)<1, object@priceOutside^(1-gamma), 0)

     sharesIn     <- meanval*prices^(1-sigma[nests])
     sharesAcross <- tapply(sharesIn,nests,sum,na.rm=TRUE)
     sharesIn     <- sharesIn / sharesAcross[nests]
     sharesAcross <- sharesAcross^((1-gamma)/(1-sigma))
     sharesAcross <- sharesAcross / (sum(sharesAcross,na.rm=TRUE) + outVal)

     shares       <- sharesIn * sharesAcross[nests]

     ##transform revenue shares to quantity shares
     if(!revenue){shares <- (shares/prices)/sum((1-sum(shares,na.rm=TRUE))/object@priceOutside,shares/prices,na.rm=TRUE)}

     names(shares) <- object@labels

     return(as.vector(shares))

}
 )



setMethod(
 f= "elast",
 signature= "CESNests",
 definition=function(object,preMerger=TRUE,market=FALSE){

     nests    <- object@nests
     gamma    <- object@slopes$gamma
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval

     shares <- calcShares(object,preMerger,revenue=TRUE)
     sharesNests <- shares/tapply(shares,nests,sum,na.rm=TRUE)[nests]

       if(market){

          alpha       <- object@slopes$alpha
          if(is.null(alpha)){
               stop("'shareInside' must be between 0 and 1 to  calculate Market Elasticity")}
          elast <- (1+alpha) * (1-gamma) * sum(shares) * (1 - sum(shares))
          names(elast) <- NULL

         }

     else{
         nprods <-  length(shares)

         elast <- diag(sigma - gamma)
         elast <- elast[nests,nests]
         elast <- elast * matrix(sharesNests,ncol=nprods,nrow=nprods,byrow=TRUE)
         elast <- elast + (gamma-1) * matrix(shares,ncol=nprods,nrow=nprods,byrow=TRUE)
         diag(elast) <- diag(elast) - sigma[nests]

         dimnames(elast) <- list(object@labels,object@labels)
     }
      return(elast)

}
 )





setMethod(
          f= "CV",
          signature= "CESNests",
          definition=function(object,revenueInside){

              alpha       <- object@slopes$alpha

                if(is.null(alpha)) stop("'shareInside' must be between 0 and 1 to  calculate Compensating Variation")

              nests       <- object@nests
              gamma       <- object@slopes$gamma
              sigma       <- object@slopes$sigma
              meanval     <- object@slopes$meanval
              shareInside <- object@shareInside


              VPre  <- sum(tapply(meanval *  object@pricePre^(1-sigma[nests]),nests,sum,na.rm=TRUE) ^((1-gamma)/(1-sigma)))
              VPost <- sum(tapply(meanval * object@pricePost^(1-sigma[nests]),nests,sum,na.rm=TRUE) ^((1-gamma)/(1-sigma)))

              ##tempPre  <- log( sum( tapply(meanval * object@pricePre^(1-sigma[nests]),nests,sum) ^((1-gamma)/(1-sigma)) ) )
              ##tempPre  <- (gamma/(gamma-1)) * log( sum( tapply( meanval^((1-gamma)/(1-sigma[nests])) * object@pricePre^(-gamma),nests,sum)^((gamma-1)/gamma)) ) - tempPre

              ##tempPost  <- log( sum( tapply(meanval * object@pricePost^(1-sigma[nests]),nests,sum) ^((1-gamma)/(1-sigma)) ) )
              ##tempPost  <- (gamma/(gamma-1)) * log( sum( tapply( meanval^((1-gamma)/(1-sigma[nests])) * object@pricePost^(-gamma),nests,sum)^((gamma-1)/gamma)) ) - tempPost


              result <- log(VPost/VPre) / ((1+alpha)*(1-gamma))
              names(result) <- NULL
              if(missing(revenueInside)){
                  warning("'revenueInside' is missing. Calculating CV as a percentage change in (aggregate) income")
                  return(result*100)}

              else{
                  totExp <- revenueInside*(1+alpha)
                  return(totExp*(exp(result)-1))
                  }

 })



ces.nests <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      nests=rep(1,length(shares)),
                      shareInside = 1,
                      normIndex=ifelse(sum(shares) < 1,NA,1),
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=1,
                      priceStart = prices,
                      isMax=FALSE,
                      constraint = TRUE,
                      parmsStart,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
                      ){



    nests <- factor(nests,levels=unique(nests))


    if(missing(parmsStart)){
        nNests <- nlevels(nests)
        parmsStart <- cumsum(runif(nNests+1,1,1.5)) # parameter values are assumed to be greater than 1

        if(constraint){parmsStart <- parmsStart[1:2]}
                            }

    ## Create CESNests  container to store relevant data
    result <- new("CESNests",prices=prices, shares=shares,margins=margins,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  nests=nests,
                  normIndex=normIndex,
                  parmsStart=parmsStart,
                  priceStart=priceStart,
                  constraint=constraint,
                  shareInside=shareInside,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)

    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,preMerger=TRUE,isMax=isMax,...)
    result@pricePost <- calcPrices(result,preMerger=FALSE,isMax=isMax,subset=subset,...)


    return(result)

}

