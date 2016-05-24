setClass(
         Class = "PCAIDS",
         contains="AIDS",
         representation=representation(

         knownElast      = "numeric",
         knownElastIndex = "numeric"
         ),

         validity=function(object){




             nprods <- length(object@shares)

             if(!(object@knownElastIndex %in% seq(1,nprods)) ){
                 stop("'knownElastIndex' value must be between 1 and the length of 'shares'")}
             if(nprods != length(object@mcDelta)){
                 stop("'mcDelta' must have the same length as 'shares'")}
             if(object@knownElast>0 || object@mktElast > 0 ){
                 stop("'mktElast', 'knownElast' must be non-positive")}
             if(abs(object@knownElast) < abs(object@mktElast) ){
                  stop("'mktElast' must be less  than 'knownElast' in absolute value")}
         }

         )








setMethod(
 f= "calcSlopes",
 signature= "PCAIDS",
 definition=function(object){


     ## Uncover linear demand slopes from shares, knownElast and mktElast
     ## Since demand is linear IN LOG PRICE, model assumes that slopes remain unchanged following merger

     shares    <- object@shares
     diversion <- object@diversion
     labels    <- object@labels

     nprod    <- length(shares)

     idx      <- object@knownElastIndex

     shareKnown <- shares[idx]


     minD <- function(betas){

       #enforce symmetry
       bknown = betas[idx]
       betas  =  betas[-idx]


        B = diag(nprod)

       B[upper.tri(B)] <- betas
       B=t(B)
       B[upper.tri(B)] <- betas
       diag(B)= 1-rowSums(B)



       m1 = bknown - shareKnown * (object@knownElast + 1 - shareKnown * (object@mktElast + 1))
       m2 = as.vector(diversion +  t(B)/diag(B)) #measure distance between observed and predicted diversion


       measure=c(m1,m2)

       return(sum(measure^2))
     }


     ## Create starting values for optimizer
     bKnown = shareKnown * (object@knownElast + 1 - shareKnown * (object@mktElast + 1))
     bStart <- bKnown*diversion[idx,]/diversion[,idx]
     bStart = -diversion*bStart
     bStart = c(bKnown,bStart[upper.tri(bStart)])


     ## create bounds for optimizer
     upper<-lower<-bStart
     upper[1]=0
     upper[-(1)]=Inf
     lower[1]=-Inf
     lower[-(1)]=0

     bestBetas=optim(bStart,minD,method="L-BFGS-B",upper=upper,lower=lower)


     B = diag(nprod)
     B[upper.tri(B)] <- bestBetas$par[-(1)]
     B=t(B)
     B[upper.tri(B)] <- bestBetas$par[-(1)]
     diag(B)= 1-rowSums(B)


     dimnames(B) <- list(labels,labels)
     object@slopes <- B
     object@intercepts <- as.vector(shares - B%*%log(object@prices))
     names(object@intercepts) <- object@labels

     return(object)
 }

      )









pcaids <- function(shares,knownElast,mktElast=-1,
                   prices,diversions,
                   ownerPre,ownerPost,
                   knownElastIndex=1,
                   mcDelta=rep(0, length(shares)),
                   subset=rep(TRUE, length(shares)),
                   priceStart=runif(length(shares)),
                   isMax=FALSE,
                   labels=paste("Prod",1:length(shares),sep=""),
                   ...){



    if(missing(prices)){ prices <- rep(NA_real_,length(shares))}

    if(missing(diversions)){
        diversions <- tcrossprod(1/(1-shares),shares)
        diag(diversions) <- -1.000000001 #correct potential floating point issue
    }

  ## Create PCAIDS container to store relevant data
    result <- new("PCAIDS",shares=shares,prices=prices,
                   quantities=shares, margins=shares,mcDelta=mcDelta,
                  subset=subset,
                  knownElast=knownElast,mktElast=mktElast,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  knownElastIndex=knownElastIndex,
                  diversion=diversions,
                  priceStart=priceStart,labels=labels)

    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate Demand Slope Coefficients
    result <- calcSlopes(result)

    ## Solve Non-Linear System for Price Changes
    result@priceDelta <- calcPriceDelta(result,isMax=isMax,subset=subset,...)


    ## Calculate Pre and Post merger equilibrium prices
    ## These are equal to NA in pcaids
    result@pricePre  <- calcPrices(result,TRUE)
    result@pricePost <- calcPrices(result,FALSE)


    return(result)
}



