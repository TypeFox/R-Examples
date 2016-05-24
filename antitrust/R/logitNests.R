
setClass(
         Class   = "LogitNests",
         contains="Logit",

         representation=representation(
         nests="factor",
         parmsStart="numeric",
         constraint="logical"
         ),

         prototype=prototype(
         parmsStart      =  numeric()
         ),

         validity=function(object){




             nprods    <- length(object@prices)
             nNestParm <- nlevels(object@nests) #calculate the number of nesting parameters

             ## Identify Singleton Nests
             nestCnt   <- tapply(object@prices,object@nests,length)
             nestCnt   <- nestCnt[object@nests]
             isSingleton <- nestCnt==1

             nNestParm <- nNestParm - sum(isSingleton) #singleton nests are not identified

             if(identical(nNestParm,1)) stop("'logit.nests', 'logit.nests.alm' may not be used for non-nested problems or problems with only singleton nests. Use 'logit', 'logit.alm' instead")

             if(nprods != length(object@nests)){
                 stop("'nests' length must equal the number of products")}



             if(!object@constraint &&
                any(tapply(object@margins[!isSingleton],object@nests[!isSingleton],
                           function(x){if(all(is.na(x))){return(TRUE)} else{return(FALSE)}}
                           )
                    ,na.rm=TRUE)
                ){
                 stop("when 'constraint' is FALSE, at least one product margin must be supplied for each non-singleton nest")
             }



             return(TRUE)
         }

         )


setMethod(
          f= "calcSlopes",
          signature= "LogitNests",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              nprods      <-   length(shares)

              margins      <-  object@margins

              ## identify which products have enough margin information
              ##  to impute Bertrand margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))

              prices       <-  object@prices
              idx          <-  object@normIndex

              parmsStart   <- object@parmsStart
              nests        <- object@nests
              nestCnt      <- tapply(prices,nests,length)
              constraint   <- object@constraint

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




              sharesNests <- tapply(shares,nests,sum)[nests]

              sharesNests <- shares / sharesNests

              revenues <- prices * shares





              ## Minimize the distance between observed and predicted margins
              minD <- function(theta){

                  alpha <- theta[1]
                  sigma <- as.numeric(isSingletonNest)
                  sigma[!isSingletonNest] <- theta[-1]

                  ## The following returns the transpose of the elasticity matrix
                  elasticity <- diag((1/sigma-1)*alpha)
                  elasticity <- elasticity[nests,nests]
                  elasticity <- elasticity * matrix(sharesNests*prices,ncol=nprods,nrow=nprods)
                  elasticity <- -1*(elasticity + alpha * matrix(shares*prices,ncol=nprods,nrow=nprods))
                  diag(elasticity) <- diag(elasticity) + (1/sigma[nests])*alpha*prices



                  elasticity <- elasticity[isMargin,isMargin]
                  revenues   <- revenues[isMargin]
                  ownerPre   <- ownerPre[isMargin,isMargin]
                  margins    <- margins[isMargin]

                  #marginsCand <- -1 * as.vector(ginv(elasticity * ownerPre) %*% (revenues * diag(ownerPre))) / revenues
                  #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

                  measure <- revenues * diag(ownerPre) + as.vector((elasticity * ownerPre) %*% (margins * revenues))
                  measure <- sum(measure^2,na.rm=TRUE)

                  return(measure)
              }

              ##  Constrained optimizer to look for solutions where alpha<0,  1 > sigma > 0.
              ##  sigma > 1 or sigma < 0 imply complements
              lowerB <- upperB <- rep(0,length(parmsStart))
              lowerB[1] <- -Inf

              upperB[-1] <- 1

              minTheta <- optim(parmsStart,minD,method="L-BFGS-B",lower= lowerB,upper=upperB)


              if(minTheta$convergence != 0){
                  warning("'calcSlopes' nonlinear solver did not successfully converge. Reason: '",minTheta$message,"'")
              }



              minAlpha           <- minTheta$par[1]
              names(minAlpha)    <- "Alpha"


              minSigma <-  as.numeric(isSingletonNest)
              minSigma[!isSingletonNest] <- minTheta$par[-1]


              minSigmaOut        <- minSigma
              minSigma           <- minSigma[nests]
              names(minSigmaOut) <- levels(nests)


                ## create index variables, contingent on whether an outside good is defined
              if(is.na(idx)){
                  idxShare <- 1 - object@shareInside
                  idxShareIn <- 1
                  idxPrice   <- object@priceOutside
                  idxSigma   <- 1

              }

              else{


                  idxShare   <- shares[idx]
                  idxShareIn <- sharesNests[idx]
                  idxPrice   <- prices[idx]
                  idxSigma   <- minSigma[idx]

               }


              meanval <-
                  log(shares) - log(idxShare) -
                      minAlpha*(prices - idxPrice) -
                          (1-minSigma)*log(sharesNests) +
                              (1-idxSigma)*log(idxShareIn)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,sigma=minSigmaOut,meanval=meanval)

              return(object)
          }
          )

setMethod(
 f= "calcShares",
 signature= "LogitNests",
 definition=function(object,preMerger=TRUE,revenue=FALSE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     nests     <- object@nests
     alpha     <- object@slopes$alpha
     sigma     <- object@slopes$sigma
     meanval   <- object@slopes$meanval
     isOutside <- as.numeric(object@shareInside < 1)
     outVal    <- ifelse(object@shareInside<1, exp(alpha*object@priceOutside), 0)

     sharesIn     <- exp((meanval+alpha*prices)/sigma[nests])

     inclusiveValue <- log(tapply(sharesIn,nests,sum,na.rm=TRUE))
     sharesAcross <-   exp(sigma*inclusiveValue)
     sharesAcross <- sharesAcross/(outVal + sum(sharesAcross,na.rm=TRUE))


     sharesIn     <- sharesIn/exp(inclusiveValue)[nests]


     shares       <- sharesIn * sharesAcross[nests]

     if(revenue){shares <- prices*shares/sum(prices*shares,object@priceOutside*(1-sum(shares,na.rm=TRUE)),na.rm=TRUE)}

     names(shares) <- object@labels

     return(as.vector(shares))

}
 )



setMethod(
 f= "elast",
 signature= "LogitNests",
 definition=function(object,preMerger=TRUE,market=FALSE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}


     nests    <- object@nests
     alpha    <- object@slopes$alpha
     sigma    <- object@slopes$sigma
     meanval  <- object@slopes$meanval

     shares <- calcShares(object,preMerger,revenue=FALSE)

     if(market){

         elast <- alpha * sum(shares*prices) * (1 - sum(shares))
         names(elast) <- NULL
         }

     else{

         sharesNests <- shares/tapply(shares,nests,sum,na.rm=TRUE)[nests]



         nprods <-  length(shares)

         elast <- diag((1/sigma-1)*alpha)
         elast <- elast[nests,nests]
         elast <- elast * matrix(sharesNests*prices,ncol=nprods,nrow=nprods,byrow=TRUE)
         elast <- -1*(elast + alpha * matrix(shares*prices,ncol=nprods,nrow=nprods,byrow=TRUE))
         diag(elast) <- diag(elast) + (1/sigma[nests])*alpha*prices

         dimnames(elast) <- list(object@labels,object@labels)

         }
      return(elast)

}
 )





setMethod(
          f= "CV",
          signature= "LogitNests",
          definition=function(object){

              nests       <- object@nests
              alpha       <- object@slopes$alpha
              sigma       <- object@slopes$sigma
              meanval     <- object@slopes$meanval



              VPre  <- sum( tapply(exp((meanval + object@pricePre*alpha)  / sigma[nests]),nests,sum,na.rm=TRUE) ^ sigma )
              VPost <- sum( tapply(exp((meanval + object@pricePost*alpha) / sigma[nests]),nests,sum,na.rm=TRUE) ^ sigma )



              result <- log(VPost/VPre)/alpha
              names(result) <- NULL
              return(result)

 })



logit.nests <- function(prices,shares,margins,
                        ownerPre,ownerPost,
                        nests=rep(1,length(shares)),
                        normIndex=ifelse(sum(shares) < 1,NA,1),
                        mcDelta=rep(0,length(prices)),
                        subset=rep(TRUE,length(prices)),
                        priceOutside=0,
                        priceStart = prices,
                        isMax=FALSE,
                        constraint = TRUE,
                        parmsStart,
                        labels=paste("Prod",1:length(prices),sep=""),
                        ...
                        ){


    nests <- factor(nests,levels = unique(nests)) # factor nests, keeping levels in the order in which they appear
    nNestParm <- sum(tapply(nests,nests,length)>1) # count the number of  non-singleton nests
    nMargins  <- length(margins[!is.na(margins)])
    maxNests  <- nMargins - 1


     if(nNestParm > maxNests){
                 stop("Additional margins must be supplied in order to calibrate nesting parameters")
             }

    if(missing(parmsStart)){

        nNests <- nlevels(nests)
        parmsStart <- runif(nNests+1) # nesting parameter values are assumed to be between 0 and 1
        parmsStart[1] <- -1* parmsStart[1] # price coefficient is assumed to be negative

        if(constraint){parmsStart <- parmsStart[1:2]}
                   }


    if(constraint && length(parmsStart)!=2){
        stop("when 'constraint' is TRUE, 'parmsStart' must be a vector of length 2")
    }
    else if(!constraint && nNestParm + 1 != length(parmsStart)){
        stop("when 'constraint' is FALSE, 'parmsStart' must be a vector of length ",nNestParm + 1)

    }


    ## Create LogitNests  container to store relevant data
    result <- new("LogitNests",prices=prices, margins=margins,
                  shares=shares,mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  nests=nests,
                  normIndex=normIndex,
                  parmsStart=parmsStart,
                  constraint=constraint,
                  priceStart=priceStart,shareInside=sum(shares),
                  labels=labels)

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

