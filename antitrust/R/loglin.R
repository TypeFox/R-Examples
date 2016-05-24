setClass(
         Class = "LogLin",
         contains="Linear",
         prototype=prototype(
         symmetry=FALSE
         ),
          validity=function(object){




             nprods <- length(object@prices)
             if(any(is.na(object@margins))){
                 stop("'margins' cannot contain NA values")
                 }
             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'prices'")}
 })


setMethod(
 f= "calcSlopes",
 signature= "LogLin",
 definition=function(object){



     margins <- object@margins
     quantities <- object@quantities
     prices <- object@prices
     ownerPre <- object@ownerPre

     revenues <- prices * quantities

     nprods <- length(margins)

     diversion <- object@diversion * tcrossprod(quantities,1/quantities)

     slopes <- matrix(margins * revenues,ncol=nprods, nrow=nprods,byrow=TRUE)

     slopes <- (revenues * diag(ownerPre)) / rowSums(slopes * diversion * ownerPre)
     slopes <- -t(slopes * diversion)

     dimnames(slopes) <- list(object@labels,object@labels)

     intercept <- as.vector(log(quantities) - slopes %*% log(prices))
     names(intercept) <-  object@labels

     object@slopes <- slopes
     object@intercepts <- intercept

     return(object)


 }
          )


setMethod(
 f= "calcPrices",
 signature= "LogLin",
 definition=function(object,preMerger=TRUE,subset,...){

     if(preMerger){
       owner <- object@ownerPre
       mc    <- object@mcPre
     }
     else{
       owner <- object@ownerPost
       mc    <- object@mcPost
     }

     nprods <- length(object@quantities)
     if(missing(subset)){
        subset <- rep(TRUE,nprods)
     }

     if(!is.logical(subset) || length(subset) != nprods ){stop("'subset' must be a logical vector the same length as 'quantities'")}


 FOC <- function(priceCand){

     if(preMerger){ object@pricePre <- priceCand}
     else{          object@pricePost <- priceCand}


     margins    <- 1 - mc/priceCand
     quantities <- calcQuantities(object,preMerger,revenue=TRUE)
     revenues   <- priceCand*quantities
     elasticities     <- t(elast(object,preMerger))

     thisFOC <- revenues * diag(owner) + as.vector((elasticities * owner) %*% (margins * revenues))
     thisFOC[!subset] <- revenues[!subset] #set quantity equal to 0 for firms not in subset


     return(thisFOC)
 }

 minResult <- BBsolve(object@priceStart,FOC,quiet=TRUE,...)

if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBSolve' reports: '",minResult$message,"'")}


 priceEst        <- minResult$par
 names(priceEst) <- object@labels
 return(priceEst)

}
 )


setMethod(
 f= "calcQuantities",
 signature= "LogLin",
 definition=function(object,preMerger=TRUE,...){


     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

     slopes    <- object@slopes
     intercept <- object@intercepts
     quantities <- exp(intercept) * apply(prices^t(slopes),2,prod)
     names(quantities) <- object@labels

     return(quantities)

}
 )




setMethod(
 f= "elast",
 signature= "LogLin",
 definition=function(object,preMerger=TRUE,market=FALSE){

    if(market){

        quantities <-  calcQuantities(object,preMerger)
        prices     <-  calcPrices(object,preMerger)
        elast      <-  sum(t(t(object@slopes * quantities) * 1/prices)) / sum(quantities) * sum(quantities * prices / sum(quantities))


       }

    else{
        elast    <- object@slopes
        dimnames(elast) <- list(object@labels,object@labels)
    }

      return(elast)

}
 )

setMethod(
 f= "CV",
 signature= "LogLin",
 definition=function(object){
     stop("CV method not currently available for 'LogLin' Class")

 })


setMethod(
 f= "calcPricesHypoMon",
 signature= "LogLin",
 definition=function(object,prodIndex){


     mc <- object@mcPre[prodIndex]
     pricePre <- object@pricePre

     calcMonopolySurplus <- function(priceCand){

         pricePre[prodIndex] <- priceCand
         object@pricePre     <- pricePre
         quantityCand        <- calcQuantities(object,TRUE)


         surplus             <- (priceCand-mc)*quantityCand[prodIndex]


         return(sum(surplus))
     }


     minResult <- optim(object@prices[prodIndex],calcMonopolySurplus,
                              method = "L-BFGS-B",lower = 0,
                              control=list(fnscale=-1))

     pricesHM <- minResult$par

     #priceDelta <- pricesHM/pricePre[prodIndex] - 1
     #names(priceDelta) <- object@labels[prodIndex]
     names(pricesHM) <- object@labels[prodIndex]

     return(pricesHM)

 })





loglinear <- function(prices,quantities,margins,diversions,
                      ownerPre,ownerPost,
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceStart=prices,
                      labels=paste("Prod",1:length(prices),sep=""),...
                     ){




    shares=quantities/sum(quantities)


    if(missing(diversions)){
        diversions <-  tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1.000000001 #correct potential floating point issue
    }


    result <- new("LogLin",prices=prices, quantities=quantities,margins=margins,
                  shares=shares,mcDelta=mcDelta, subset=subset, priceStart=priceStart,
                  ownerPre=ownerPre,diversion=diversions,
                  ownerPost=ownerPost, labels=labels)


    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)


    ## Calculate pre and post merger equilibrium prices
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,subset=subset,...)


   return(result)

}

