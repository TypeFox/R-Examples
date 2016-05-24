setClass(
         Class   = "LogitALM",
         contains="Logit",
         representation=representation(
          parmsStart="numeric"
         ),
         prototype=prototype(
         normIndex         =  1
         ),

         validity=function(object){



             nMargins  <- length(object@margins[!is.na(object@margins)])

             if(nMargins<2){stop("At least 2 elements of 'margins' must not be NA in order to calibrate demand parameters")}

             if(object@shareInside!=1){
                 stop(" sum of 'shares' must equal 1")
             }

              if(length(object@parmsStart)!=2){
                 stop("'parmsStart' must a vector of length 2")
                 }
         }
         )

setMethod(
          f= "calcSlopes",
          signature= "LogitALM",
          definition=function(object){

              ## Uncover Demand Coefficents

              ownerPre     <-  object@ownerPre
              shares       <-  object@shares
              margins      <-  object@margins
              prices       <-  object@prices

              nprods <- length(object@shares)

              ##identify which products have enough margin information
              ##  to impute Bertrand margins
              isMargin    <- matrix(margins,nrow=nprods,ncol=nprods,byrow=TRUE)
              isMargin[ownerPre==0]=0
              isMargin    <- !is.na(rowSums(isMargin))

              minD <- function(theta){

                  alpha <- theta[1]
                  sOut  <- theta[2]

                  probs <- shares * (1 - sOut)
                  elast <- -alpha *  matrix(prices * probs,ncol=nprods,nrow=nprods)
                  diag(elast) <- alpha*prices - diag(elast)

                  revenues <- probs * prices
                  #marginsCand <- -1 * as.vector(ginv(elast * ownerPre) %*% (revenues * diag(ownerPre))) / revenues
                  #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

                  elast      <-   elast[isMargin,isMargin]
                  revenues   <-   revenues[isMargin]
                  ownerPre   <-   ownerPre[isMargin,isMargin]
                  margins    <-   margins[isMargin]

                  #marginsCand <- -1 * as.vector(ginv(elasticity * ownerPre) %*% (revenues * diag(ownerPre))) / revenues
                  #measure <- sum((margins - marginsCand)^2,na.rm=TRUE)

                  measure <- revenues * diag(ownerPre) + as.vector((elast * ownerPre) %*% (margins * revenues))
                  measure <- sum(measure^2,na.rm=TRUE)

                  return(measure)
              }

               ## Constrain optimizer to look  alpha <0,  0 < sOut < 1
              lowerB <- c(-Inf,0)
              upperB <- c(-1e-10,.99999)

              minTheta <- optim(object@parmsStart,minD,method="L-BFGS-B",lower= lowerB,upper=upperB)$par

              if(isTRUE(all.equal(minTheta[[2]],0))){stop("Estimated outside share is close to 0. Use `logit' function instead")}
              
              meanval <- log(shares * (1 - minTheta[2])) - log(minTheta[2]) - minTheta[1] * (prices - object@priceOutside)

              names(meanval)   <- object@labels


              object@slopes      <- list(alpha=minTheta[1],meanval=meanval)
              object@shareInside <- 1-minTheta[2]

              return(object)

          }

          )


logit.alm <- function(prices,shares,margins,
                      ownerPre,ownerPost,
                      mcDelta=rep(0,length(prices)),
                      subset=rep(TRUE,length(prices)),
                      priceOutside=0,
                      priceStart = prices,
                      isMax=FALSE,
                      parmsStart,
                      labels=paste("Prod",1:length(prices),sep=""),
                      ...
                      ){


    if(missing(parmsStart)){
        parmsStart <- runif(2)
        parmsStart[1] <- -1* parmsStart[1] # price coefficient is assumed to be negative
    }

    ## Create Logit  container to store relevant data
    result <- new("LogitALM",prices=prices, shares=shares,
                  margins=margins,
                  ownerPre=ownerPre,
                  ownerPost=ownerPost,
                  mcDelta=mcDelta,
                  subset=subset,
                  priceOutside=priceOutside,
                  priceStart=priceStart,
                  shareInside=sum(shares),
                  parmsStart=parmsStart,
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


