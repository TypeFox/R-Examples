
setClass(
         Class   = "LogitCap",
         contains="Logit",
         representation=representation(
         mktSize               = "numeric",
         capacities           = "numeric"

         ),


         validity=function(object){





             nprods <- length(object@shares)


             if(nprods != length(object@capacities)){
                 stop("'prices', 'capacities' must all be vectors with the same length")}

             if(any(is.na(object@capacities) |
                    !is.finite(object@capacities) |
                    object@capacities<0 ,na.rm=TRUE)){stop("'capacities' values must be positive, finite numbers")}

             if(length(object@mktSize)!=1 || isTRUE(object@mktSize<0)){stop("mktSize must be a positive number")}

             if(any(object@mktSize*object@shares > object@capacities)){stop("utilization is greater than capacity")}

             if(identical(object@mktSize*object@shares,object@capacities)){stop("utilization cannot equal capacity for all products")}

             if(any(is.na(object@margins[object@mktSize*object@shares == object@capacities]))){
                 stop("'margins' cannot equal NA for capacity constrained products")
                 }

                return(TRUE)

         })


setMethod(
          f= "calcSlopes",
          signature= "LogitCap",
          definition=function(object){

              ## Uncover Demand Coefficents


              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices
              capacities  <-  object@capacities/object@mktSize
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

              ##create a matrix of 1s and 0s where the i,jth element equals 1 if product i is NOT producing at capacity
              notBinds <- matrix(as.numeric(capacities > shares),ncol=nprods,nrow=nprods,byrow=TRUE)
              ## create a TRUE/FALSE vector equal to TRUE if a single product firm is capacity constrained
              singleConstrained <- rowSums( object@ownerPre>0) == 1 & capacities == shares

              ## Minimize the distance between observed and predicted margins
              minD <- function(alpha){

                  ## the following returns the elasticity TRANSPOSED
                  elast <- -alpha *  matrix(prices * shares,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices + diag(elast)


                  FOC <- revenues * diag(ownerPre) + (elast * ownerPre * notBinds) %*% (margins * revenues)

                  ## omit the FOCs of single product, capacity constrained firms
                  measure <- sum(as.vector(FOC[!singleConstrained])^2,na.rm=TRUE)

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
 f= "calcQuantities",
 signature= "LogitCap",
 definition=function(object,preMerger=TRUE){

     quantities <- object@mktSize * calcShares(object,preMerger,revenue=FALSE)

     return(quantities)

 }
          )

## compute margins
setMethod(
 f= "calcMargins",
 signature= "LogitCap",
 definition=function(object,preMerger=TRUE){


     if( preMerger) {
         margins <- object@margins #capacity-constrained margins not identified -- use supplied margins
         constrained <- object@capacities == object@mktSize*object@shares

         owner  <- object@ownerPre
         revenue<- calcShares(object,preMerger,revenue=TRUE)[!constrained]
         elast <-  elast(object,preMerger)
         margins[!constrained] <-  -1 * as.vector(ginv(t(elast*owner)[!constrained,!constrained]) %*% revenue) / revenue

     }

     else{
         prices <- object@pricePost
         mc     <- object@mcPost

         margins <- 1 - mc/prices
     }


     names(margins) <- object@labels

     return(as.vector(margins))
     }

 )



setMethod(
 f= "calcPrices",
 signature= "LogitCap",
 definition=function(object,preMerger=TRUE,isMax=FALSE,subset,...){


     capacities <- object@capacities

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
         capacities <- capacities[subset]
         }


     priceEst <- rep(NA,nprods)

     ##Define system of FOC as a function of prices
     FOC <- function(priceCand){

         if(preMerger){ object@pricePre[subset] <- priceCand}
         else{          object@pricePost[subset] <- priceCand}


         margins          <- 1 - mc/priceCand
         quantities       <- calcQuantities(object,preMerger)[subset]
         revenues         <- quantities * priceCand
         elasticities     <- elast(object,preMerger)[subset,subset]

         thisFOC <- revenues * diag(owner) + as.vector(t(elasticities * owner) %*% (margins * revenues))
         constraint <- quantities - capacities

         measure <- thisFOC + constraint + sqrt(thisFOC^2 + constraint^2)

         return(measure)
     }


     ## Find price changes that set FOCs equal to 0
     minResult <- BBsolve(object@priceStart,FOC,quiet=TRUE,...)

      if(minResult$convergence != 0){warning("'calcPrices' nonlinear solver may not have successfully converged. 'BBsolve' reports: '",minResult$message,"'")}

     priceEst[subset]        <- minResult$par
     names(priceEst) <- object@labels

  if(isMax){

         hess <- genD(FOC,minResult$par) #compute the numerical approximation of the FOC hessian at optimium
         hess <- hess$D[,1:hess$p]
         hess <- hess * (owner>0)   #0 terms not under the control of a common owner

         state <- ifelse(preMerger,"Pre-merger","Post-merger")

         if(any(eigen(hess)$values>0)){warning("Hessian of first-order conditions is not positive definite. ",state," price vector may not maximize profits. Consider rerunning 'calcPrices' using different starting values")}
     }

     return(priceEst)

 }
 )



setMethod(
 f= "calcPricesHypoMon",
 signature= "LogitCap",
 definition=function(object,prodIndex,...){


     mc       <- object@mcPre[prodIndex]
     pricePre <- object@pricePre

      FOC <- function(priceCand){

          thisPrice <- pricePre
          thisPrice[prodIndex] <- priceCand

          object@pricePre <- thisPrice

          margins          <- 1 - mc/priceCand
          quantities       <- calcQuantities(object,preMerger=TRUE)[prodIndex]
          revenues         <- quantities * priceCand
          elasticities     <- elast(object,preMerger=TRUE)[prodIndex,prodIndex]

          thisFOC <- revenues + as.vector(t(elasticities) %*% (margins * revenues))
          constraint <- quantities - object@capacities[prodIndex]

          measure <- thisFOC + constraint + sqrt(thisFOC^2 + constraint^2)

         return(measure)
      }



      ## Find price changes that set FOCs equal to 0
     minResult <- BBsolve(object@priceStart[prodIndex],FOC,quiet=TRUE,...)

     if(minResult$convergence != 0){warning("'calcPricesHypoMon' nonlinear solver may not have successfully converged. 'BBSolve' reports: '",minResult$message,"'")}


     pricesHM <- minResult$par
      #priceDelta <- pricesHM/pricePre[prodIndex] - 1
      #names(priceDelta) <- object@labels[prodIndex]
     names(priceHM) <- object@labels[prodIndex]

     return(priceHM)

 })



logit.cap <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      capacities,
                      mktSize=sum(capacities),
                      normIndex=ifelse(sum(shares)<1,NA,1),
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=0,
                      priceStart = prices,
                      isMax=FALSE,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
                      ){


    ## Create LogitCap  container to store relevant data
    result <- new("LogitCap",prices=prices, shares=shares,
                  margins=margins,capacities=capacities, mktSize=mktSize,
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

