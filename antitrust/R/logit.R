
setClass(
         Class   = "Logit",
         contains="Bertrand",
         representation=representation(

         prices           = "numeric",
         margins          = "numeric",
         priceStart       = "numeric",
         normIndex        = "vector",
         shareInside      = "numeric",
         priceOutside     = "numeric"

         ),
    prototype=prototype(
      priceStart  = numeric(),
      normIndex   = 1,
      shareInside = numeric(),
      priceOutside = 0
    ),
         validity=function(object){



             margins <- object@margins

             nprods <- length(object@shares)
             sumShares <- sum(object@shares)

             if(
                 nprods != length(margins) ||
                 nprods != length(object@prices)){
                 stop("'prices', 'margins' and 'shares' must all be vectors with the same length")}

             if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")


             if(any(margins<0 | margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")

             if(!(is.matrix(object@ownerPre))){
                 ownerPre <- ownerToMatrix(object,TRUE)
             }
             else{ownerPre <- object@ownerPre}


              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))

             if(!any(isMargin)) stop("Insufficient margin information to calibrate demand parameters.")

             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'shares'")}

              if(
                !(object@shareInside >=0 &&
                  object@shareInside <=1)
                ){
                 stop("'shareInside' must be between 0 and 1")
         }

              if(
                 !(all(object@shares >0) &&
                   all(object@shares <=1))
                 ){
                 stop("elements of vector 'shares' must be between 0 and 1")
         }

             if(!(
                  (isTRUE(all.equal(sumShares,1))  && object@normIndex %in% 1:nprods) ||
                   (sumShares < 1 && is.na(object@normIndex))
                  )){
                 stop("'normIndex' must take on a value between 1 and ",nprods,
                      " if 'shares' sum to 1 , or NA if the sum of shares is less than 1")
             }

             if(length(object@priceOutside) != 1 || object@priceOutside<0
                ){stop("'priceOutside' must be a non-negative number")}

             return(TRUE)

         })



setMethod(
          f= "calcSlopes",
          signature= "Logit",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex

              if(is.na(idx)){
                  idxShare <- 1 - object@shareInside
                  idxPrice <- object@priceOutside
              }
              else{
                  idxShare <- shares[idx]
                  idxPrice <- prices[idx]
               }

              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)
              revenues <- shares * prices

              ## identify which products have enough margin information
              ##  to impute Bertrand margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))


              ## Minimize the distance between observed and predicted margins
              minD <- function(alpha){

                  ## the following returns the elasticity TRANSPOSED
                  elast <- -alpha *  matrix(prices * shares,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices + diag(elast)


                  elast      <- elast[isMargin,isMargin]
                  revenues   <- revenues[isMargin]
                  ownerPre   <- ownerPre[isMargin,isMargin]
                  margins    <- margins[isMargin]

                  #marginsCand <- -1 * as.vector(ginv(elast * ownerPre) %*% (revenues * diag(ownerPre))) / revenues
                  #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

                  measure <- revenues * diag(ownerPre) + as.vector((elast * ownerPre) %*% (margins * revenues))
                  measure <- sum(measure^2,na.rm=TRUE)

                  return(measure)
              }

              minAlpha <- optimize(minD,c(-1e6,0))$minimum


              meanval <- log(shares) - log(idxShare) - minAlpha * (prices - idxPrice)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=minAlpha,meanval=meanval)


              return(object)
          }
          )



setMethod(
 f= "calcPrices",
 signature= "Logit",
 definition=function(object,preMerger=TRUE,isMax=FALSE,subset,...){


     priceStart <- object@priceStart

     if(preMerger){
       owner <- object@ownerPre
       mc    <- object@mcPre
     }
     else{
       owner <- object@ownerPost
       mc    <- object@mcPost
     }

     nprods <- length(object@shares)
     if(missing(subset)){
        subset <- rep(TRUE,nprods)
     }

     if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'shares'")}

     if(any(!subset)){
         owner <- owner[subset,subset]
         mc    <- mc[subset]
         priceStart <- priceStart[subset]
         }

     priceEst <- rep(NA,nprods)



     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         if(preMerger){ object@pricePre[subset] <- priceCand}
         else{          object@pricePost[subset] <- priceCand}


         margins   <- 1 - mc/priceCand
         revenues  <- calcShares(object,preMerger,revenue=TRUE)[subset]
         elasticities     <- t(elast(object,preMerger)[subset,subset])

         thisFOC <- revenues * diag(owner) + as.vector((elasticities * owner) %*% (margins * revenues))

         return(thisFOC)
     }

     ## Find price changes that set FOCs equal to 0
     minResult <- BBsolve(priceStart,FOC,quiet=TRUE,...)

      if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}


     if(isMax){

         hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
         hess <- hess$D[,1:hess$p]
         hess <- hess * (owner>0)   #0 terms not under the control of a common owner

         state <- ifelse(preMerger,"Pre-merger","Post-merger")

         if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. ",state," price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
     }


     priceEst[subset]        <- minResult$par
     names(priceEst) <- object@labels

     return(priceEst)

 }
 )


setMethod(
 f= "calcShares",
 signature= "Logit",
 definition=function(object,preMerger=TRUE,revenue=FALSE){


     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}


     alpha    <- object@slopes$alpha
     meanval  <- object@slopes$meanval

     outVal <- ifelse(object@shareInside<1, exp(alpha*object@priceOutside), 0)

     shares <- exp(meanval + alpha*prices)
     shares <- shares/(outVal+ sum(shares,na.rm=TRUE))

     if(revenue){shares <- prices*shares/sum(prices*shares,object@priceOutside*(1-sum(shares,na.rm=TRUE)),na.rm=TRUE)}

     names(shares) <- object@labels

     return(shares)

}
 )



setMethod(
 f= "elast",
 signature= "Logit",
 definition=function(object,preMerger=TRUE,market=FALSE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}


     labels <- object@labels

     alpha    <- object@slopes$alpha

     shares <-  calcShares(object,preMerger)

     if(market){

         elast <- alpha * sum(shares*prices,na.rm=TRUE) * (1 - sum(shares,na.rm=TRUE))

         }

     else{



         nprods <-  length(shares)

         elast <- -alpha  * matrix(prices*shares,ncol=nprods,nrow=nprods,byrow=TRUE)
         diag(elast) <- alpha*prices + diag(elast)

         dimnames(elast) <- list(labels,labels)
     }

      return(elast)

 }
          )



setMethod(
          f= "CV",
          signature= "Logit",
          definition=function(object){

              alpha       <- object@slopes$alpha
              meanval     <- object@slopes$meanval


              VPre  <- sum(exp(meanval + object@pricePre*alpha),na.rm=TRUE)
              VPost <- sum(exp(meanval + object@pricePost*alpha),na.rm=TRUE)

              result <- log(VPost/VPre)/alpha

              return(result)

 })



setMethod(
 f= "calcPricesHypoMon",
 signature= "Logit",
 definition=function(object,prodIndex){


     mc       <- object@mcPre[prodIndex]
     pricePre <- object@pricePre

     calcMonopolySurplus <- function(priceCand){

         pricePre[prodIndex] <- priceCand #keep prices of products not included in HM fixed at premerger levels
         object@pricePre     <- pricePre
         sharesCand          <- calcShares(object,TRUE,revenue=FALSE)

         surplus             <- (priceCand-mc)*sharesCand[prodIndex]

         return(sum(surplus,na.rm=TRUE))
     }


     maxResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                              method = "L-BFGS-B",lower = 0,
                              control=list(fnscale=-1))

     pricesHM <- maxResult$par

     #priceDelta <- pricesHM/pricePre[prodIndex] - 1
     #names(priceDelta) <- object@labels[prodIndex]
     names(pricesHM) <- object@labels[prodIndex]

     return(pricesHM)

 })




logit <- function(prices,shares,margins,
                  ownerPre,ownerPost,
                  normIndex=ifelse(sum(shares)<1,NA,1),
                  mcDelta=rep(0,length(prices)),
                  subset=rep(TRUE,length(prices)),
                  priceOutside = 0,
                  priceStart = prices,
                  isMax=FALSE,
                  labels=paste("Prod",1:length(prices),sep=""),
                  ...
                  ){


    ## Create Logit  container to store relevant data
    result <- new("Logit",prices=prices, shares=shares,
                  margins=margins,
                  normIndex=normIndex,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
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

