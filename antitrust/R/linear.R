setClass(
         Class = "Linear",
         contains="Bertrand",
         representation=representation(
         intercepts       = "vector",
         prices           = "vector",
         quantities       = "numeric",
         margins          = "numeric",
         priceStart       = "numeric",
         diversion        = "matrix",
         symmetry         = "logical"
         ),
          prototype=prototype(
          intercepts    =  numeric(),
          symmetry      =  TRUE
        ),
         validity=function(object){



             nprods <- length(object@shares) # count the number of products
             diversion <- object@diversion

             if(nprods != length(object@quantities) ||
                nprods != length(object@margins) ||
                nprods != length(object@prices)){
                 stop("'prices', 'quantities', 'margins', and 'shares' must all be vectors with the same length")}

             if(any(object@prices<0,na.rm=TRUE))             stop("'prices' values must be positive")
             if(any(object@quantities<0,na.rm=TRUE))          stop("'quantities' values must be positive")
             if(any(object@margins<0 | object@margins>1,na.rm=TRUE)) stop("'margins' values must be between 0 and 1")


             if(!isTRUE(all.equal(diag(diversion),rep(-1,nprods)))){ stop("'diversions' diagonal elements must all equal -1")}

             diag(diversion)=1
             if(any(diversion > 1 | diversion<0)){
                 stop("'diversions' off-diagonal elements must be between 0 and 1")}

            
             if (any(rowSums(object@diversion,na.rm=TRUE)>0,na.rm=TRUE)){ stop("'diversions' rows cannot sum to greater than 0")}

             if(nprods != nrow(object@diversion) ||
                nprods != ncol(object@diversion)){
                 stop("'diversions' must be a square matrix")
             }


             if(nprods != length(object@priceStart)){
                 stop("'priceStart' must have the same length as 'shares'")}


             if(!is.logical(object@symmetry) || length(object@symmetry)!=1){stop("'symmetry' must equal TRUE or FALSE")}

             if(!object@symmetry &&
                length(object@margins[!is.na(object@margins)])!= nprods){
                 stop("When 'symmetry' is FALSE, all product margins must be supplied")
                 }

             return(TRUE)

         }
 )


setMethod(
 f= "calcSlopes",
 signature= "Linear",
 definition=function(object){

     margins    <- object@margins
     quantities <- object@quantities
     prices     <- object@prices
     diversion  <- object@diversion
     ownerPre   <- object@ownerPre
     symmetry  <- object@symmetry

     nprod <- length(margins)



      if(!symmetry){


          slopes <- matrix(margins * prices,ncol=nprod, nrow=nprod,byrow=TRUE)
          slopes <- diag(ownerPre)/rowSums(slopes * diversion * ownerPre) * quantities
          slopes <- -t(slopes * diversion)


      }

     else{

       existMargins <- which(!is.na(margins))

       revenues <- prices*quantities
       k <- existMargins[1] ## choose a diagonal demand parameter corresponding to a provided margin


       minD <- function(betas){

         #enforce symmetry

         B=diag(betas[1:nprod])
         B[upper.tri(B,diag=FALSE)] <- betas[-(1:nprod)]
         B=t(B)
         B[upper.tri(B,diag=FALSE)] <- betas[-(1:nprod)]
      
         elast <- t(B * tcrossprod(1/quantities,prices))

         marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% (revenues * diag(ownerPre))) / revenues


         m1 <- margins - marginsCand
         m2 <- as.vector(diversion +  t(B)/diag(B)) #measure distance between observed and predicted diversion


         measure=c(m1,m2)

         return(sum(measure^2,na.rm=TRUE))
       }


       ## Create starting values for optimizer

       bKnown      =  -quantities[k]/(prices[k]*margins[k])
       bStart      =   bKnown*diversion[k,]/diversion[,k]
       bStart      =  -diversion*bStart
       parmStart   =   c(diag(bStart),bStart[upper.tri(bStart,diag=FALSE)])



       ## constrain diagonal elements so that D'b >=0
       ## constrain off-diagonal elements to be non-negative.
       
       ui          =  diag(length(parmStart))
       ui[1:nprod,1:nprod] = t(diversion)
       
       ci = rep(0,length(parmStart)) 
       
       
       
       bestParms=constrOptim(parmStart,minD,grad=NULL,ui=ui,ci=ci)

       slopes = diag(bestParms$par[1:nprod])

       slopes[upper.tri(slopes,diag=FALSE)] <- bestParms$par[-(1:nprod)]
       slopes=t(slopes)
       slopes[upper.tri(slopes,diag=FALSE)] <- bestParms$par[-(1:nprod)]
      

}


     dimnames(slopes) <- list(object@labels,object@labels)


     intercept <- as.vector(quantities - slopes %*% prices)
     names(intercept) <- object@labels

     if(!symmetry &&
        !isTRUE(all.equal(slopes,t(slopes)))){
         warning("Matrix of demand slopes coefficients is not symmetric. Demand parameters may not be consistent with utility maximization theory.")}

     if(any(intercept<0))   {warning(  "Some demand intercepts are negative")}
     if(any(diag(slopes)>0)){warning(  "Some own-slope coefficients are positive")}

     object@slopes <- slopes
     object@intercepts <- intercept

     return(object)


 }
          )




setMethod(
 f= "calcPrices",
 signature= "Linear",
 definition=function(object,preMerger=TRUE,subset,...){

     slopes    <- object@slopes
     intercept <- object@intercepts
     priceStart<- object@priceStart

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



     ##first try the analytic solution
#      prices <-
#          solve((slopes*diag(owner)) + (t(slopes)*owner)) %*%
#              ((t(slopes)*owner) %*% mc - (intercept*diag(owner)))
#
#       prices <- as.vector(prices)
#       quantities <- as.vector(intercept + t(slopes) %*% prices)
#
#      ##use the numeric solution if analytic solution yields negative quantities
#      if(any(subset) || any(quantities<0)){
#
#          if(any(quantities<0)){
#           warning("Equilibrium prices yield negative equilibrium quantities. Recomputing equilibrium prices  under the restriction that equilbrium quantities must be non-negative")
#                               }
#          else if(any(subset)){
#            warning("Elements of 'subset' are  FALSE. Computing equilbrium under the restriction that these products have 0 sales")
#          }

         FOC <- function(priceCand){

             if(preMerger){ object@pricePre  <- priceCand}
             else{          object@pricePost <- priceCand}

             margins   <- priceCand - mc
             quantities  <- calcQuantities(object,preMerger)

             thisFOC <- quantities*diag(owner) + (t(slopes)*owner) %*% margins
             thisFOC[!subset] <- quantities[!subset] #set quantity equal to 0 for firms not in subset

             return(as.vector(crossprod(thisFOC)))

         }

         ##Find starting value that always meets boundary conditions
         ##startParm <- as.vector(solve(slopes) %*% (-intercept + 1))

         minResult <- constrOptim(object@priceStart,FOC,grad=NULL,ui=slopes,ci=-intercept,...)
         
         if(!isTRUE(all.equal(minResult$convergence,0))){
             warning("'calcPrices' solver may not have successfully converged.'constrOptim' reports: '",minResult$message,"'")
           }

         prices <- minResult$par

     #}

     names(prices) <- object@labels

     return(prices)


 }
          )



setMethod(
 f= "calcQuantities",
 signature= "Linear",
 definition=function(object,preMerger=TRUE){

     if(preMerger){ prices <- object@pricePre}
     else{          prices <- object@pricePost}

      slopes    <- object@slopes
      intercept <- object@intercepts


     quantities <- as.vector(intercept+slopes %*% prices)
     names(quantities) <- object@labels

     return(quantities)

}
 )

setMethod(
 f= "calcShares",
 signature= "Linear",
 definition=function(object,preMerger=TRUE,revenue=FALSE){

     quantities <- calcQuantities(object,preMerger)

     if (revenue){
         if(preMerger){ prices <- object@pricePre}
         else{          prices <- object@pricePost}

         return(prices*quantities/sum(prices*quantities))
     }

     else{return(quantities/sum(quantities))}
 }
          )


setMethod(
 f= "elast",
 signature= "Linear",
 definition=function(object,preMerger=TRUE,market=FALSE){

       if(preMerger){ prices <- object@pricePre}
       else{          prices <- object@pricePost}

       slopes    <- object@slopes


       quantities <-  calcQuantities(object,preMerger)

       if(market){

           elast <-sum(slopes)/sum(quantities) * sum(quantities * prices / sum(quantities))

       }

       else{


           elast <- slopes * tcrossprod(1/quantities,prices)
           dimnames(elast) <- list(object@labels,object@labels)
       }

       return(elast)

}
 )

setMethod(
 f= "CV",
 signature= "Linear",
 definition=function(object){

     slopes    <- object@slopes

     if(!isTRUE(all.equal(slopes,t(slopes)))){
                  stop("price coefficient matrix must be symmetric in order to calculate compensating variation. Suggest setting 'symmetry=TRUE'")
              }

     intercept <- object@intercepts
     pricePre  <- object@pricePre
     pricePost <- object@pricePost

     result <- sum(intercept*(pricePost-pricePre)) + .5 * as.vector(pricePost%*%slopes%*%pricePost - pricePre%*%slopes%*%pricePre)

     return(result)
 })


setMethod(
 f= "calcPricesHypoMon",
 signature= "Linear",
 definition=function(object,prodIndex){

     nprods <- length(prodIndex)
     intercept <- object@intercepts
     slopes <- object@slopes
     mc <- object@mcPre[prodIndex]
     pricePre <- object@pricePre

     calcMonopolySurplus <- function(priceCand){


         pricePre[prodIndex] <- priceCand
         quantityCand <- intercept + as.vector(slopes %*% pricePre)

         surplus <- (priceCand-mc)*quantityCand[prodIndex]

         return(sum(surplus))
     }

     ##Find starting value that always meets boundary conditions
     ##Note: if nprods=1, need to use a more accurate optimizer.

     if(nprods > 1){

         if(det(slopes)!=0){startParm <- as.vector(solve(slopes) %*% (1 - intercept ))}
         else{startParm <- rep(0,nprods)}


         priceConstr <- pricePre
         priceConstr[prodIndex] <- 0

         maxResult <- constrOptim(startParm[prodIndex],calcMonopolySurplus,
                                  grad=NULL,
                                  ui=slopes[prodIndex,prodIndex],
                                  ci=-intercept[prodIndex] - as.vector(slopes %*% priceConstr)[prodIndex],
                                  control=list(fnscale=-1))

         pricesHM <- maxResult$par
     }


     else{

         upperB <- -(intercept[prodIndex] + sum(pricePre[-prodIndex]*slopes[prodIndex,-prodIndex]))/slopes[prodIndex,prodIndex]

         maxResult <- optimize(calcMonopolySurplus,c(0,upperB),maximum = TRUE)
         pricesHM <- maxResult$maximum
      }

     #priceDelta <- pricesHM/pricePre[prodIndex] - 1
     #names(priceDelta) <- object@labels[prodIndex]
     names(pricesHM) <- object@labels[prodIndex]

     return(pricesHM)


 })



linear <- function(prices,quantities,margins, diversions, symmetry=TRUE,
                   ownerPre,ownerPost,
                   mcDelta=rep(0,length(prices)),
                   subset=rep(TRUE,length(prices)),
                   priceStart=prices,
                   labels=paste("Prod",1:length(prices),sep=""),
                   ...
                     ){

    shares <- quantities/sum(quantities)

    if(missing(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1.000000001 #correct potential floating point issue
        
    }


     result <- new("Linear",prices=prices, quantities=quantities,margins=margins,
                   shares=shares,mcDelta=mcDelta, subset=subset,
                   ownerPre=ownerPre,diversion=diversions, symmetry=symmetry,
                   ownerPost=ownerPost, priceStart=priceStart,labels=labels)


     ## Convert ownership vectors to ownership matrices
     result@ownerPre  <- ownerToMatrix(result,TRUE)
     result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients and Intercepts
    result <- calcSlopes(result)


    ## Calculate marginal cost
    result@mcPre <-  calcMC(result,TRUE)
    result@mcPost <- calcMC(result,FALSE)

    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,subset=subset,...)


   return(result)

}


