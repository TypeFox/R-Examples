sim <- function(prices,demand=c("Linear","AIDS","LogLin","Logit","CES","LogitNests","CESNests","LogitCap"),demand.param,
                ownerPre,ownerPost,nests, capacities,
                mcDelta=rep(0,length(prices)),
                subset=rep(TRUE,length(prices)),
                priceOutside,
                priceStart,
                labels=paste("Prod",1:length(prices),sep=""),...){

    demand <- match.arg(demand)
    nprods <- length(prices)

    if(missing(priceStart)){
        if(demand=="AIDS"){priceStart <- runif(nprods)}
        else{              priceStart <- prices}
        }

    ## Create placeholders values to fill required Class slots

    shares <- margins <- rep(1/nprods,nprods)


    if(!missing(nests)){nests <- factor(nests,levels=unique(nests))}


    ## general checks
    if(!is.list(demand.param)){stop("'demand.param' must be a list.")}

    ## Checks for discrete choice models
    if(demand %in% c("CESNests","LogitNests","CES","Logit","LogitCap")){

         if(!("meanval" %in% names(demand.param))){
                stop("'demand.param' does not contain 'meanval'.")
            }
         if(length(demand.param$meanval) != nprods || any(is.na(demand.param$meanval))){
             stop("'meanval' must be a length-k vector of product mean valuations. NAs not allowed.")
             }

         if(demand %in% c("LogitNests","Logit","LogitCap")){

             ## An outside option is assumed to exist if all mean valuations are non-zero
             if(all(demand.param$meanval!=0)){
                 normIndex <- NA
                 shares <- rep(1/(nprods+1),nprods)
             }
             else{
                 normIndex <- which(demand.param$meanval==0)

                 if(length(normIndex)>1){
                     warning("multiple values of meanval are equal to zero. Normalizing with respect to the first product with zero mean value.")
                     normIndex <- normIndex[1]
                 }

             }

             if(!("alpha" %in% names(demand.param))   ||
                  length(demand.param$alpha) != 1     ||
                  isTRUE(demand.param$alpha>0)){
                 stop("'demand.param' does not contain 'alpha' or 'alpha' is not a negative number.")
             }

             shareInside <- sum(shares)
             if(missing(priceOutside)){priceOutside <- 0}

         }


         else  if(demand %in% c("CESNests","CES")){
             if(!("gamma" %in% names(demand.param))   ||
                  length(demand.param$gamma) != 1     ||
                  isTRUE(demand.param$gamma<0)){
                 stop("'demand.param' does not contain 'gamma' or 'gamma' is not a positive number.")
             }


             ## uncover Numeraire Coefficients
             if(!("alpha" %in% names(demand.param)) &&
                !("shareInside" %in% names(demand.param))){
                 warning("'demand.param' does not contain either 'alpha' or 'shareInside'. Setting shareInside=1 and alpha=NULL.")
                 shareInside=1
                 demand.param$alpha=NULL
             }

             else if("shareInside" %in% names(demand.param)){
                 shareInside=demand.param$shareInside
                 demand.param$shareInside <- NULL

                 if(shareInside<1) {demand.param$alpha <- 1/shareInside -1}
                 else{ demand.param$alpha <- NULL}


             }


             ## An outside option is assumed to exist if all mean valuations are non-zero
             if(all(demand.param$meanval!=1)){
                 normIndex <- NA
                 shares <- rep(1/(nprods+1),nprods)
             }
             else{
                 normIndex <- which(demand.param$meanval==1)

                 if(length(normIndex)>1){
                     warning("multiple values of meanval are equal to one. Normalizing with respect to the first product with  mean value equal to 1.")
                     normIndex <- normIndex[1]
                 }

             }


             if(missing(priceOutside)){priceOutside <- 1}
         }

         if(demand %in% c("CESNests","LogitNests")){

             if(!("sigma" %in% names(demand.param))){
                 stop("'demand.param' does not contain 'sigma'.")
            }
             if(length(demand.param$sigma)==1){constraint=TRUE}
             else{constraint=FALSE}


             if(missing(nests) ||
                length(nests)!= nprods ){stop("When 'demand' equals 'CESNests' or 'LogitNests', 'nests' must equal a vector whose length equals the number of products.")}

             if(nlevels(nests) != length(demand.param$sigma) &&
                length(demand.param$sigma) != 1){
                 stop("The number of nests in 'nests' must either equal 1 or the number of nesting parameters in 'demand.param$sigma'.")}

         }


     }


    ## Checks for Linear-demand style models
    if(demand %in% c("Linear","LogLin","AIDS")){

        if(!("slopes" %in% names(demand.param))){stop("'demand.param' does not contain 'slopes'")}
        if(!("intercepts" %in% names(demand.param))){stop("'demand.param' does not contain 'intercepts'")}

        if(!(is.matrix(demand.param$slopes))   ||
           ncol(demand.param$slopes)!=nprods   ||
           nrow(demand.param$slopes)!=nprods   ||
           any(diag(demand.param$slopes)>0)){
            stop("'slopes' must be a k x k matrix of slope coeficients whose diagonal elements must all be negative.")}
        if(!is.vector(demand.param$intercepts)     ||
           length(demand.param$intercepts)!=nprods ||
           isTRUE(any(demand.param$intercepts<0,na.rm=TRUE))){
            stop("'intercepts' must be a length-k vector whose elements are all non-negative")
        }

        if (demand == "AIDS" &&
            !("mktElast" %in% names(demand.param))){
            warning("'demand.param' does not contain 'mktElast'. Setting 'mktElast' equal to -1")
            demand.param$mktElast=-1

        }

    }






    ## Create constructors for each demand system specified in the 'demand' parameter

    if(demand == "CESNests"){

         result <- new(demand,prices=prices, shares=shares,margins=margins,
                       mcDelta=mcDelta,
                       subset=subset,
                       ownerPre=ownerPre,
                       ownerPost=ownerPost,
                       nests=nests,
                       normIndex=normIndex,
                       parmsStart=c(demand.param$gamma,demand.param$sigma),
                       priceStart=priceStart,
                       constraint=constraint,
                       shareInside=shareInside,labels=labels)

    }

    else if(demand == "LogitNests"){

        result <- new(demand,prices=prices, shares=shares,margins=margins,
                      mcDelta=mcDelta,
                      subset=subset,
                      ownerPre=ownerPre,
                      ownerPost=ownerPost,
                      nests=nests,
                      normIndex=normIndex,
                      parmsStart=c(demand.param$alpha,demand.param$sigma),
                      priceStart=priceStart,
                      constraint=constraint,
                      shareInside=shareInside,labels=labels)

    }


    else if(demand %in% c("Logit","CES")){


        result <- new(demand,prices=prices, shares=shares,
                      margins=margins,
                      normIndex=normIndex,
                      mcDelta=mcDelta,
                      subset=subset,
                      ownerPre=ownerPre,
                      ownerPost=ownerPost,
                      priceStart=priceStart,shareInside=shareInside,
                      labels=labels)

    }


    else if(demand == "LogitCap"){

        if(!("mktSize" %in% names(demand.param))){
            if(!missing(capacities) ){
                warning("'demand.param' does not contain 'mktSize'. Setting 'mktSize' equal to the sum of 'capacities'.")
                mktSize <- sum(capacities)
            }
            else{stop("'demand.param' does not contain 'mktSize'")}
        }
        else{mktSize <- demand.param$mktSize}


        shares <- capacities/mktSize
        shares <- shares/sum(shares)

        result <- new(demand, prices=prices, shares=shares,
                      margins=margins,capacities=capacities, mktSize=mktSize,
                      normIndex=normIndex,
                      ownerPre=ownerPre,
                      ownerPost=ownerPost,
                      mcDelta=mcDelta,
                      subset=subset,
                      priceStart=priceStart,shareInside=shareInside,
                      labels=labels)
    }


    else if(demand == "Linear"){



        result <- new(demand,prices=prices, quantities=shares,margins=margins,
                      shares=shares,mcDelta=mcDelta,  subset=subset,
                      ownerPre=ownerPre,diversion=-diag(nprods),
                      symmetry=identical(demand.param$slopes,t(demand.param$slopes)),
                      ownerPost=ownerPost, priceStart=priceStart,labels=labels)

    }

    else if(demand == "AIDS"){

        ## find the market elasticity that best explains user-supplied intercepts and prices

        aidsShares    <- as.vector(demand.param$intercepts + demand.param$slopes %*% log(prices)) # AIDS needs actual shares for prediction
        aidsDiv       <- tcrossprod(1/(1-aidsShares),aidsShares)
        diag(aidsDiv) <- -1

        result <- new(demand,prices=prices, quantities=shares,margins=margins,
                      shares=aidsShares,
                      mcDelta=mcDelta,  subset=subset,mktElast=demand.param$mktElast,
                      ownerPre=ownerPre,diversion=aidsDiv,
                      priceStart=priceStart,
                      ownerPost=ownerPost, labels=labels)

    }



    else if(demand == "LogLin"){


        result <- new(demand,prices=prices, quantities=shares,margins=margins,
                      shares=shares,mcDelta=mcDelta, subset=subset, priceStart=priceStart,
                      ownerPre=ownerPre,diversion=-diag(nprods),
                      ownerPost=ownerPost, labels=labels)

    }


    if(demand %in% c("Linear","LogLin","AIDS")){
        result@slopes <- demand.param$slopes
        result@intercepts <- demand.param$intercepts
    }
    else{result@slopes=demand.param}


    ## Convert ownership vectors to ownership matrices
    result@ownerPre  <- ownerToMatrix(result,TRUE)
    result@ownerPost <- ownerToMatrix(result,FALSE)

    ## Calculate marginal cost
    result@mcPre     <-  calcMC(result,TRUE)
    result@mcPost    <-  calcMC(result,FALSE)

    if(demand == "AIDS"){
        ## Solve Non-Linear System for Price Changes
        result@priceDelta <- calcPriceDelta(result,...)
    }


    ## Solve Non-Linear System for Price Changes
    result@pricePre  <- calcPrices(result,TRUE,...)
    result@pricePost <- calcPrices(result,FALSE,subset=subset,...)


    return(result)
}
