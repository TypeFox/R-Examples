##Consider switching the order in which parameters are solved for: put distribution parms first


setClass(
         Class   = "Auction2ndCap",
         contains="Antitrust",
         representation=representation(
         capacities       = "numeric",
         margins          = "numeric",
         prices           = "numeric",
         reserve          = "numeric",
         shareInside      = "numeric",
         sellerCostCDF    = "character",
         sellerCostCDFLowerTail    = "logical",
         sellerCostPDF    = "function",
         sellerCostBounds = "numeric",
         sellerCostParms  = "numeric",
	       buyerValuation        = "numeric", #was buyerValuation
         reservePre       = "numeric",
         reservePost      = "numeric",
         mcDelta          = "numeric",
         parmsStart       = "numeric"

         ),
         prototype=prototype(
         reservePre      =  numeric(),
         reserveost     =  numeric(),
         buyerValuation  =  numeric(),
         sellerCostParms =  numeric(),
         sellerCostCDFLowerTail    = TRUE


          ),
         validity=function(object){

             nprods <- length(object@labels)

             cdf    <- object@sellerCostCDF
             if(is.na(object@reserve)){parmsStart <- object@parmsStart[-1]}
             else{parmsStart <- object@parmsStart}

             if(nprods != length(object@capacities) ||
                nprods != length(object@margins)    ||
                nprods != length(object@prices)     ||
                nprods != length(object@ownerPre)   ||
                nprods != length(object@ownerPost)  ||
                nprods != length(object@mcDelta)
                ){
                 stop("'capacities', 'margins', 'prices', 'ownerPre', 'ownerPost', 'mcDelta' must all be the same length")}
             if(length(object@reserve) != 1 ||
                (!is.na(object@reserve) && object@reserve<0)){
                 stop("'reserve'must be a length 1 vector whose value is greater than 0")}
             if(
                 (!is.na(object@shareInside) &&
                  (object@shareInside > 1 ||
                   object@shareInside < 0)) ||
                 length(object@shareInside) != 1){
                 stop("'shareInside' must be a length 1 vector whose value is between 0 and 1")
             }
             if(any(object@capacities<0 | is.na(object@capacities),na.rm=TRUE)){
                 stop("'capacities' cannot be negative or equal to NA")
             }
             if(any(object@prices<=0,na.rm=TRUE)){
                 stop("'prices' must be positive")}
             if(any(object@margins > 1 | object@margins < 0,na.rm=TRUE)){
                 stop("'margins' must be between 0 and 1")
             }


             if( identical(cdf,"punif")){
                 if(length(parmsStart)!=2){
                     if(is.na(object@reserve)){stop("parmsStart must be a length 3 vector whose first element is the starting value for 'reserve'")}
                     else{stop("For the Uniform distribution, 'parmsStart' must be a numeric vector of length 2")}
                 }
                    if(parmsStart[1] >= parmsStart[2]){
                  stop("The upper bound must be greater than the lower bound")
                }
             }

             else if( identical(cdf,"pexp")){
                 if(length(parmsStart)!=1){
                     if(is.na(object@reserve)){stop("parmsStart must be a length 2 vector whose first element is the starting value for 'reserve'")}
                     else{stop("For the Exponential distribution, 'parmsStart' must be a numeric vector of length 1")}
                 }
                 if(parmsStart[1] <=0){
                    stop("For the Exponential distribution, 'parmsStart' must be a numeric vector of length ",
                         1 + is.na(object@reserve),"whose final element is greater than 0")
                 }
             }

             else if( identical(cdf,"pweibull")){
                   if(length(parmsStart)!=2){
                       if(is.na(object@reserve)){stop("parmsStart must be a length 3 vector whose first element is the starting value for 'reserve'")}
                       else{stop("For the Weibull distribution, 'parmsStart' must be a numeric vector of length 2")}
                     }
                   if(parmsStart[1] <=0  ||
                      parmsStart[2] <=0 ){
                       stop("For the Weibull distribution, 'parmsStart' must be a numeric vector of length ",
                            2 + is.na(object@reserve)," whose final 2 elements must be greater than 0")
                   }
               }

             else if( identical(cdf,"pgumbel")){
                 if(length(parmsStart)!=2){
                     if(is.na(object@reserve)){stop("parmsStart must be a length 3 vector whose first element is the starting value for 'reserve'")}
                     else{stop("For the Gumbel distribution, 'parmsStart' must be a numeric vector of length 2")}
                     }
                 if(parmsStart[2] <=0){
                     stop("For the Gumbel distribution, 'parmsStart' must be a numeric vector of length ",
                          2 + is.na(object@reserve)," whose final element must be greater than 0")
                 }
             }
             else if( identical(cdf,"pfrechet")){
                if(length(parmsStart)!=3){
                    if(is.na(object@reserve)){stop("parmsStart must be a length 4 vector whose first element is the starting value for 'reserve'")}
                    else{stop("For the Frechet distribution, 'parmsStart' must be a numeric vector of length 3")}
                }
                if(parmsStart[2] <=0  ||
                   parmsStart[3] < 2  ){
                      stop("For the Frechet distribution, 'parmsStart' must be a numeric vector of length ",
                           2 + is.na(object@reserve)," whose next-to-last element must be positive and whose final element must be at least 2")
                  }
            }


             data  <- sum(!is.na(object@reserve) , !is.na(object@shareInside) , sum(!is.na(object@prices)) , sum(!is.na(object@margins)))
             unknowns <-  is.na(object@reserve) + length(parmsStart)

                if(data < unknowns ){
                    stop("Insufficient information to calibrate model parameters: ",unknowns, " unknowns but only ", data, " pieces of information")
                }

             return(TRUE)
         }
         )

setGeneric (
 name= "calcOptimalReserve",
 def=function(object,...){standardGeneric("calcOptimalReserve")}
 )
setGeneric (
 name= "cdfG",
 def=function(object,...){standardGeneric("cdfG")}
 )
setGeneric (
  name= "calcExpectedPrice",
  def=function(object,...){standardGeneric("calcExpectedPrice")}
)
setGeneric (
 name= "calcExpectedLowestCost",
 def=function(object,...){standardGeneric("calcExpectedLowestCost")}
 )
setGeneric (
 name= "calcBuyerExpectedCost",
 def=function(object,...){standardGeneric("calcBuyerExpectedCost")}
 )
setGeneric (
 name= "calcProducerSurplus",
 def=function(object,...){standardGeneric("calcProducerSurplus")}
 )
setGeneric (
  name= "calcSellerCostParms",
  def=function(object,...){standardGeneric("calcSellerCostParms")}
)
setGeneric (
  name= "calcBuyerValuation",
  def=function(object,...){standardGeneric("calcBuyerValuation")}
)

setMethod(
    f= "calcSellerCostParms",
    signature= "Auction2ndCap",
    definition=function(object,...){

        ## method to calibrate seller cost distribution parameters
        sellerCostParms <- object@sellerCostParms
        reserve         <- object@reserve
        shareInside     <- object@shareInside
        margins         <- object@margins
        prices          <- object@prices
        parmsStart      <- object@parmsStart


        cdf <- object@sellerCostCDF

        minD <- function(parmsStart){

            if(is.na(reserve)){r <- parmsStart[1]; parmsStart <- parmsStart[-1]}
            else{r <- reserve}

            sellerCostParms <- parmsStart

            object@reservePre      <- r
            object@sellerCostParms <- sellerCostParms

            ## For uniform, frechet,  distribution bounds are a function
            ## of distribution parameters

            if(identical(cdf,"punif")){
                object@sellerCostBounds <- parmsStart }
            else if(identical(cdf,"pfrechet")){
               object@sellerCostBounds[1] <- parmsStart[1] }





            ##calculate each bidder's profit margin, conditional on bidder winning
            thisInShare <- cdfG(object,preMerger=TRUE)
            thisMargin  <- calcProducerSurplus(object,preMerger=TRUE,exAnte=FALSE)
            thisPrice   <- calcPrices(object,preMerger=TRUE,exAnte=FALSE)


            measure   <- c(margins - thisMargin/thisPrice,
                           1 - thisPrice/prices,
                           shareInside - thisInShare)

            return(sum(measure^2,na.rm=TRUE))

        }


        if(identical(cdf,"punif")) {

         ui = diag(length(parmsStart))
         if(is.na(reserve)){ui[1,nrow(ui)-1]=-1} #constrain reserve to be greater than cLower
         ui[nrow(ui),nrow(ui)-1]=-1 #constrain cLower to be less than cUpper
         ci = rep(0,length(parmsStart))


         result <- constrOptim(parmsStart,minD,grad=NULL,ui=ui,ci=ci,...)

       }
        else{
          lb <- ub <- rep(Inf,length(parmsStart))

          if( identical(cdf,"pexp")){lb[1] <- 1e-20}
          else if( identical(cdf,"pweibull")){ lb[1:2] <- 1e-20} #shape,scale must be positive
          else if( identical(cdf,"pgumbel")){  lb[2] <- 1e-20} #scale must be positive
          else if( identical(cdf,"pfrechet")){ lb[2] <- 1e-20; lb[3] <- 2}  #scale must be positive, shape must be > 2 for finite variance

          if(is.na(reserve)){lb <- c(0,lb); ub <- c(Inf,ub)}

          if(length(parmsStart)>1){method="L-BFGS-B"}
          else{method="Brent"; ub=1e12} #'Brent' is equivalent to using optimize for 1D problems

          result <- optim(parmsStart,minD,method=method,lower=lb,upper=ub,...)

        }

        result <- result$par
        if(is.na(reserve)){object@reserve <- result[1]; result <- result[-1]}
        object@sellerCostParms <- result

        if(identical(cdf, "punif")){
          object@sellerCostBounds <- result }
        else if(identical(cdf,"pfrechet")){
          object@sellerCostBounds[1] <- result[1] }


        return(object)
  }
)

setMethod(
  f= "calcBuyerValuation",
  signature= "Auction2ndCap",
  definition=function(object){

    ## Use FOC from buyers cost minimization problem
    ## to uncover buyer cost parameter
    capacities <- object@capacities
    totCap       <- sum(capacities)
    reserve    <- object@reserve

    object@reservePre <- reserve
    cdfF = match.fun(object@sellerCostCDF)
    pdfF = object@sellerCostPDF

    sellerCostParms <- c(list(reserve),as.list(object@sellerCostParms))
    fr = do.call(pdfF,sellerCostParms)

    sellerCostParms <- c(sellerCostParms,
                         lower.tail=as.list(object@sellerCostCDFLowerTail))
    Fr = do.call(cdfF,sellerCostParms)

    gr <- totCap*fr*(1-Fr)^(totCap-1)

    expectedPrice  <- calcExpectedPrice(object,preMerger=TRUE)
    partialSupplierProfits <- (1-Fr)^(totCap-capacities) - (1-Fr)^totCap
    partialSupplierProfits <- sum(partialSupplierProfits)/gr

    result <- reserve + partialSupplierProfits

    return(result)

  }
)
setMethod(
          f= "calcOptimalReserve",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE,lower,upper){

            if(missing(lower)){lower <- max(object@sellerCostBounds[1],0)}
            if(missing(upper)){upper <- object@buyerValuation}

            minD <- function(r){
              if(preMerger){object@reservePre <- r}
              else{object@reservePost <- r}
                  calcBuyerExpectedCost(object,preMerger=preMerger)
            }

            res <- optimize(
              f  = minD,
              lower = lower,
              upper = upper,
            )


              rStar <- res$minimum
              return(rStar)
          }
)


setMethod(
          f= "calcBuyerExpectedCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){


              shareInside <- cdfG(object,preMerger=preMerger)
              val  <- object@buyerValuation * (1-shareInside) + calcExpectedPrice(object,preMerger=preMerger)*shareInside

              return(val)
          })

setMethod(
          f= "calcPrices",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE,exAnte=TRUE){


            val <- calcProducerSurplus(object,preMerger=preMerger,exAnte=exAnte) + calcMC(object,preMerger=preMerger,exAnte=exAnte)

            names(val) <- object@labels
            return(val)
          })

setMethod(
          f= "cdfG",
          signature= "Auction2ndCap",
          definition=function(object,c,preMerger=TRUE){

              if(missing(c)){
                  if(preMerger){
                      c <- object@reservePre
                      capacities <- sum(object@capacities)
                            }
                  else{
                      c <- object@reservePost
                      capacities <- sum(object@capacities*(1+object@mcDelta))
                   }


              }

              else{

                  if(preMerger){capacities <- object@capacities}
                  else{capacities <- tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum)}


              }


              cdfF = match.fun(object@sellerCostCDF)
              sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                                   lower.tail=as.list(object@sellerCostCDFLowerTail))

              Fc = do.call(cdfF,sellerCostParms)
              retval = 1-(1-Fc)^capacities

              if(!preMerger && length(capacities)>1){

                  temp <- rep(NA, length(object@ownerPre))
                  temp[object@ownerPre == object@ownerPost] <- retval
                  retval <- temp

              }
                  return(retval)

              }
          )

setMethod(
  f= "calcExpectedPrice",
  signature= "Auction2ndCap",
  definition=function(object,preMerger=TRUE){

    val <- calcExpectedLowestCost(object,preMerger=preMerger) + sum(calcProducerSurplus(object,preMerger=preMerger),na.rm=TRUE)/cdfG(object,preMerger=preMerger)
    return(val)
  }
)

setMethod(
  f= "calcMC",
  signature= "Auction2ndCap",
  definition=function(object,t,preMerger=TRUE,exAnte=TRUE){


    cdfF <- match.fun(object@sellerCostCDF)
    pdfF <- object@sellerCostPDF
    sellerCostBounds <-object@sellerCostBounds



    if(preMerger) {
      capacities <- object@capacities
      r    <- object@reservePre
    }
    else {
      capacities <- tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum)
      r    <- object@reservePost
    }

    totCap <- sum(capacities)

    if(missing(t)){t <- capacities}



    ## The expected production cost
    ecIntegrand = function(c,t){
      sellerCostParms <- c(list(c),as.list(object@sellerCostParms))

      fc = do.call(pdfF,sellerCostParms)

      sellerCostParms <- c(sellerCostParms,
                           lower.tail=as.list(object@sellerCostCDFLowerTail))
      Fc = do.call(cdfF,sellerCostParms)

      retval = t*c*fc*(1-Fc)^(totCap-1)
      retval = ifelse(is.finite(retval),retval,0)

      return(retval)
    }

    result <- sapply(
      t,
      function(t.i) {
        if( r < sellerCostBounds[2]) {
            retval = integrate(ecIntegrand,lower=sellerCostBounds[1],upper=r, stop.on.error = FALSE,t=t.i)$value
        }
        else {
          retval = integrate(ecIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2], stop.on.error = FALSE,t=t.i)$value
        }

        return(retval)
      })


    if(!preMerger && length(t)>1){

        temp <- rep(NA, length(object@ownerPre))
        temp[object@ownerPre == object@ownerPost] <- result
        result <- temp

        }

    if(!exAnte){result <- result/calcShares(object,preMerger=preMerger,exAnte=TRUE)}

    return(result)

  }
)


setMethod(
          f= "calcExpectedLowestCost",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE){


              if(preMerger){capacities <- object@capacities}
              else{capacities <- object@capacities*(1+object@mcDelta)}

              num    <- calcMC(object,t=sum(capacities),preMerger=preMerger,exAnte=TRUE)


              retval <- num/cdfG(object, preMerger=preMerger)
              return(retval)


          })


setMethod(
          f= "calcProducerSurplus",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE,exAnte=TRUE){

             sellerCostBounds <-object@sellerCostBounds
              if(preMerger){r    <- object@reservePre}
              else{r    <- object@reservePost}


              if(preMerger) { capacities = object@capacities }
              else {          capacities = tapply(object@capacities*(1+object@mcDelta),object@ownerPost,sum) }

              totCap = sum(capacities)

              espIntegrand = function(c,t){
                  sellerCostParms <- c(list(c),as.list(object@sellerCostParms),
                                       lower.tail=as.list(object@sellerCostCDFLowerTail))
                  Fc <- do.call(match.fun(object@sellerCostCDF),sellerCostParms)
                  val <- (1-Fc)^(totCap-t)-(1-Fc)^totCap
              }


              retval <- sapply(
                                capacities,
                                function(t.i) {

                                  if( r < sellerCostBounds[2]) {
                                      retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=r,
                                                          stop.on.error = FALSE,t=t.i)$value
                                  }
                                  else {
                                    retval <- integrate(espIntegrand,lower=sellerCostBounds[1],upper=sellerCostBounds[2],
                                                        stop.on.error = FALSE,t=t.i)$value
                                  }

                                    return(retval)
                                })

             if(!preMerger){

                 temp <- rep(NA, length(object@ownerPre))
                 temp[object@ownerPre == object@ownerPost] <- retval
                 retval <- temp

             }

             if(!exAnte){retval <- retval/calcShares(object,preMerger=preMerger,exAnte=TRUE)}

              return(retval)
          })

setMethod(
          f= "calcShares",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE,exAnte=TRUE){


              ownerPre  <- object@ownerPre
              ownerPost  <- object@ownerPost

              shareInside <- cdfG(object,preMerger=preMerger)

              if(preMerger){
                  capacities <- object@capacities/sum(object@capacities)
                  names(capacities) <- object@labels
                  if(exAnte){return((shareInside*capacities))}
                  else{return(capacities)}
              }

              else{

                  capacities <- object@capacities*(1+object@mcDelta)
                  result <- rep(NA, length(object@ownerPre))
                  names(result) <- object@labels
                  result[object@ownerPre == ownerPost] <- tapply(capacities,ownerPost,sum)

                  if(exAnte){return((shareInside*result))}
                  else{return(result)}

              }


          }
)

setMethod(
          f= "calcMargins",
          signature= "Auction2ndCap",
          definition=function(object,preMerger=TRUE,exAnte=TRUE){

              result <- calcProducerSurplus(object,preMerger=preMerger,exAnte=exAnte)/calcPrices(object,preMerger=preMerger,exAnte=exAnte)
              return(result)
          })


##summarize method

setMethod(
 f= "summary",
 signature= "Auction2ndCap",
 definition=function(object,exAnte=FALSE,parameters=FALSE,digits=2){

     curWidth <-  getOption("width")


     pricePre   <-  calcPrices(object,preMerger=TRUE,exAnte=exAnte)
     pricePost  <-  calcPrices(object,preMerger=FALSE,exAnte=exAnte)
     priceDelta <- (pricePost/pricePre - 1) * 100


     outPre  <-  calcShares(object,TRUE,exAnte=exAnte) * 100
     outPost <-  calcShares(object,FALSE,exAnte=exAnte) * 100





     mcDelta <- object@mcDelta

     outDelta <- (outPost/outPre - 1) * 100


     isParty <- object@ownerPost != object@ownerPre
     isParty <- c(object@ownerPre[isParty],object@ownerPost[isParty])
     isParty <- factor(ifelse(object@ownerPre %in% isParty,1,0),levels=0:1,labels=c(" ","*"))

     results <- data.frame(pricePre=pricePre,pricePost=pricePost,
                           priceDelta=priceDelta,sharesPre=outPre,
                           sharesPost=outPost,sharesDelta=outDelta)


     if(sum(abs(mcDelta))>0) results <- cbind(results,mcDelta=mcDelta)


     rownames(results) <- paste(isParty,object@labels)


     cat("\nMerger simulation results under '",class(object),"':\n\n",sep="")

     options("width"=100) # this width ensures that everything gets printed on the same line
     print(round(results,digits),digits=digits)
     options("width"=curWidth) #restore to current width

     cat("\n\tNotes: '*' indicates merging parties. Deltas are percent changes.\n")

     if(exAnte){cat("\tEx Ante shares and prices are reported.\n")}
     else{cat("\tShares and prices conditional on a firm winning are reported.\n")}

     results <- cbind(isParty, results)

     cat("\n\nPre-Merger Buyer Reserve:",round(object@reservePre,digits),sep="\t")
     cat("\nPost-Merger Buyer Reserve:",round(object@reservePost,digits),sep="\t")
     cat("\n\n% Change In Expected Price:",round((calcExpectedPrice(object,FALSE)-calcExpectedPrice(object,TRUE))/calcExpectedPrice(object,TRUE)*100,digits),sep="\t")
     cat("\n")
     cat("% Change In Buyer's Expected Cost:",round((calcBuyerExpectedCost(object,FALSE)-calcBuyerExpectedCost(object,TRUE))/calcBuyerExpectedCost(object,TRUE)*100,digits),sep="\t")
     cat("\n\n")


     if(parameters){

         cat("\nSupplier Cost Distribution Parameters:\n\n")

         print(round(object@sellerCostParms,digits))

         cat("\nBuyer Valuation:\n\n")
         print(round(object@buyerValuation,digits))

             cat("\n\n")

         }

     rownames(results) <- object@labels
     return(invisible(results))

 })




## Create Constructor Function

auction2nd.cap <- function(capacities, margins,prices,reserve=NA,shareInside=NA,
                           sellerCostCDF=c("punif","pexp","pweibull","pgumbel","pfrechet"),
                           ownerPre,ownerPost,
                           mcDelta=rep(0,length(capacities)),
                           constrain.reserve=TRUE, parmsStart,
                           labels=as.character(ownerPre),...
                          ){


    sellerCostCDF <- match.arg(sellerCostCDF)
    lower.tail    <- TRUE
    sellerCostPDF <- match.fun(paste("d",substring(sellerCostCDF,2,),sep=""))

    if(is.na(reserve)){reserve <- NA_real_}
    if(is.na(shareInside)){shareInside <- NA_real_}

    if (missing(parmsStart)){
      avgPrice =  mean(prices,na.rm=TRUE)
      minPrice =  min(prices,na.rm=TRUE)
      maxPrice =  max(prices,na.rm=TRUE)


      if(identical(sellerCostCDF,"punif")){
          if(maxPrice > minPrice){
              parmsStart <- sellerCostBounds <- c(minPrice,maxPrice)} #uniform on price range
          else{parmsStart <- sellerCostBounds <- c(1e-20,maxPrice)}
      }
      else if(identical(sellerCostCDF,"pexp")){parmsStart=1/avgPrice; }
      else if(identical(sellerCostCDF,"pweibull")){ parmsStart=c(avgPrice,avgPrice)}
      else if(identical(sellerCostCDF,"pgumbel")){parmsStart=c(avgPrice,sqrt(avgPrice)) }
      else if(identical(sellerCostCDF,"pfrechet")){parmsStart=c(0,sqrt(avgPrice),2) }


      if(is.na(reserve) || missing(reserve)){
        parmsStart<-c(maxPrice,parmsStart)
      }
    }

    ##remove any names from parmStart
    names(parmsStart) <- NULL

    if(identical(sellerCostCDF,"punif")){sellerCostBounds <- parmsStart[-1]}
    else if(identical(sellerCostCDF,"pexp")){sellerCostBounds=c(0,Inf)}
    else if(identical(sellerCostCDF,"pweibull")){ sellerCostBounds=c(0,Inf)}
    else if(identical(sellerCostCDF,"pgumbel")){ sellerCostBounds=c(-Inf,Inf) ;lower.tail=FALSE}
    else if(identical(sellerCostCDF,"pfrechet")){ sellerCostBounds=c(parmsStart[2],Inf) ; lower.tail=FALSE}



    result <- new("Auction2ndCap",capacities=capacities,
                  margins=margins,
                  prices=prices,
                  reserve=reserve,
                  shareInside=shareInside,
                  sellerCostCDF=sellerCostCDF,
                  sellerCostCDFLowerTail=lower.tail,
                  sellerCostBounds=sellerCostBounds,
                  sellerCostPDF=sellerCostPDF,
                  ownerPre=ownerPre,ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  parmsStart=parmsStart,
                  labels=labels)


    ## Calibrate seller cost parameters
    result                 <- calcSellerCostParms(result,...)

    ## Calibrate buyer cost parameter
    result@buyerValuation       <- calcBuyerValuation(result)

    ## Compute pre- and post-merger reserves
    result@reservePre      <- calcOptimalReserve(result,preMerger=TRUE) #Find Buyer Reserve pre-merger
    if(constrain.reserve){ result@reservePost <- result@reservePre}
    else{result@reservePost     <- calcOptimalReserve(result,preMerger=FALSE)} #Find Buyer Reserve post-merger

    ## Compute equilbrium prices
    result@pricePre         <- calcPrices(result,preMerger=TRUE,exAnte=FALSE)
    result@pricePost        <- calcPrices(result,preMerger=FALSE,exAnte=FALSE)

    ## Compute equilibrium marginal costs
    result@mcPre            <- calcMC(result,preMerger=TRUE,exAnte=FALSE)
    result@mcPost           <- calcMC(result,preMerger=FALSE,exAnte=FALSE)

    return(result)
    }



