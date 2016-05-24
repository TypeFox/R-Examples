setClass(
    Class   = "CES",
    contains="Logit",
    prototype=prototype(
    priceOutside=0
    )
    )


setMethod(
          f= "calcSlopes",
          signature= "CES",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              idx          <-  object@normIndex
              shareInside  <-  object@shareInside

              ## uncover Numeraire Coefficients
              if(shareInside <= 1 && shareInside>0) {alpha <- 1/shareInside - 1}
              else{alpha <- NULL}

              ## if sum of shares is less than 1, add numeraire
               if(is.na(idx)){
                  idxShare <- 1 - sum(shares)
                  idxPrice <- object@priceOutside
              }
              else{
                  idxShare <- shares[idx]
                  idxPrice <- prices[idx]
               }


              ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              ## identify which products have enough margin information
              ##  to impute Bertrand margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))





              ## Minimize the distance between observed and predicted margins
              minD <- function(gamma){


                  elasticity <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods)
                  diag(elasticity) <- -gamma + diag(elasticity)

                  elasticity <- elasticity[isMargin,isMargin]
                  shares     <- shares[isMargin]
                  ownerPre   <- ownerPre[isMargin,isMargin]
                  margins    <- margins[isMargin]

                  #marginsCand <- -1 * as.vector(ginv(elasticity * ownerPre) %*% (shares * diag(ownerPre))) / shares
                  #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)
                   FOC <- (shares * diag(ownerPre)) + (elasticity * ownerPre) %*% (shares * margins)
                   measure<-sum(FOC^2,na.rm=TRUE)

                  return(measure)
              }

              minGamma <- optimize(minD,c(1,1e6))$minimum


              meanval <- log(shares) - log(idxShare) + (minGamma - 1) * (log(prices) - log(idxPrice))
              meanval <- exp(meanval)

              names(meanval)   <- object@labels

              object@slopes    <- list(alpha=alpha,gamma=minGamma,meanval=meanval)


              return(object)
          }
          )


setMethod(
 f= "calcShares",
 signature= "CES",
 definition=function(object,preMerger=TRUE,revenue=FALSE){




     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}


     gamma    <- object@slopes$gamma
     meanval  <- object@slopes$meanval

     outVal <- ifelse(sum(object@shares)<1, object@priceOutside^(1-gamma), 0)
     shares <- meanval*prices^(1-gamma)
     shares <- shares/(sum(shares,na.rm=TRUE) + outVal)

     ##transform revenue shares to quantity shares
     if(!revenue){shares <- (shares/prices)/sum((1-sum(shares,na.rm=TRUE))/object@priceOutside,shares/prices,na.rm=TRUE)}

     names(shares) <- object@labels

     return(as.vector(shares))

}
 )





setMethod(
 f= "elast",
 signature= "CES",
 definition=function(object,preMerger=TRUE,market=FALSE){

     gamma    <- object@slopes$gamma

     shares <-  calcShares(object,preMerger,revenue=TRUE)


      if(market){

          alpha       <- object@slopes$alpha
          if(is.null(alpha)){
              stop("'shareInside' must be between 0 and 1 to  calculate Market Elasticity")}
          elast <- (1+alpha) * (1-gamma) * sum(shares) * (1 - sum(shares))

         }

     else{

         nprods <-  length(shares)
         elast <- (gamma - 1 ) * matrix(shares,ncol=nprods,nrow=nprods,byrow=TRUE)
         diag(elast) <- -gamma + diag(elast)

         dimnames(elast) <- list(object@labels,object@labels)
     }
      return(elast)

}
 )



setMethod(
          f= "CV",
          signature= "CES",
          definition=function(object,revenueInside){

              alpha       <- object@slopes$alpha


             if(is.null(alpha)) stop("'shareInside' must be between 0 and 1 to  calculate Compensating Variation")

              gamma       <- object@slopes$gamma
              meanval     <- object@slopes$meanval
              shareInside <- object@shareInside



              VPre  <- sum(meanval * object@pricePre^(1-gamma),na.rm=TRUE)
              VPost <- sum(meanval * object@pricePost^(1-gamma),na.rm=TRUE)

              result <- log(VPost/VPre) / ((1+alpha)*(1-gamma))

              if(missing(revenueInside)){
                  warning("'revenueInside' is missing. Calculating CV as a percentage change in (aggregate) income")
                  return(result*100)}

              else{
                  totExp <- revenueInside*(1+alpha)
                  return(totExp*(exp(result)-1))
              }


 })


ces <- function(prices,shares,margins,
                ownerPre,ownerPost,
                shareInside = 1,
                normIndex=ifelse(sum(shares)<1,NA,1),
                mcDelta=rep(0,length(prices)),
                subset=rep(TRUE,length(prices)),
                priceOutside=1,
                priceStart = prices,
                isMax=FALSE,
                labels=paste("Prod",1:length(prices),sep=""),
                ...
                ){




    ## Create CES  container to store relevant data
    result <- new("CES",prices=prices, shares=shares, margins=margins,
                  normIndex=normIndex,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  priceStart=priceStart,
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

