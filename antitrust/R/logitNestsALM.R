setClass(
         Class   = "LogitNestsALM",
         contains="LogitNests",
         prototype=prototype(
         normIndex         =  1
         ),
          validity=function(object){


              if(!identical(object@shareInside,1)){
                 stop(" sum of 'shares' must equal 1")
             }


         }
         )

setMethod(
          f= "calcSlopes",
          signature= "LogitNestsALM",
          definition=function(object){

              ## Uncover Demand Coefficents

              ownerPre     <-  object@ownerPre
              shares       <-  object@shares

              margins      <-  object@margins
              prices       <-  object@prices

              nprods       <-  length(shares)

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
                  parmsStart   <- parmsStart[c(TRUE,TRUE,!isSingletonNest)] #always retain first two elements; these are
                                                                          # the initial value for price coefficient, outside share
              }

               ## Uncover price coefficient and mean valuation from margins and revenue shares


              nprods <- length(shares)

              sharesNests <- tapply(shares,nests,sum)[nests]

              sharesNests <- shares / sharesNests



              minD <- function(theta){

                  alpha <- theta[1]
                  sOut  <- theta[2]
                  sigma <- as.numeric(isSingletonNest)
                  sigma[!isSingletonNest] <- theta[-c(1,2)]

                  probs <- shares * (1 - sOut)
                  elast <- diag((1/sigma-1)*alpha)
                  elast <- elast[nests,nests]
                  elast <- elast * matrix(sharesNests*prices,ncol=nprods,nrow=nprods)
                  elast <- -1*(elast + alpha * matrix(probs*prices,ncol=nprods,nrow=nprods))
                  diag(elast) <- diag(elast) + (1/sigma[nests])*alpha*prices


                  revenues <- probs * prices
                  marginsCand <- -1 * as.vector(solve(elast * ownerPre) %*% (revenues * diag(ownerPre))) / revenues

                  measure <- marginsCand-margins
                  measure <- sum((measure)^2,na.rm=TRUE)

                  return(measure)
              }

               ## Constrain optimizer to look  alpha <0, 0 < sOut < 1, sigma

              lowerB <- upperB <- rep(0,length(parmsStart))
              lowerB[1] <- -Inf
              upperB[-1] <- 1

              minTheta <- optim(parmsStart,minD,method="L-BFGS-B",lower= lowerB,upper=upperB)


              minAlpha           <- minTheta$par[1]
              names(minAlpha)    <- "Alpha"

              shareOut           <- minTheta$par[2]


              minSigma <-  as.numeric(isSingletonNest)
              minSigma[!isSingletonNest] <- minTheta$par[-c(1,2)]


              minSigmaOut        <- minSigma
              minSigma           <- minSigma[nests]
              names(minSigmaOut) <- levels(nests)


               meanval <-
                  log(shares * (1 - shareOut)) - log(shareOut) -
                      minAlpha*(prices - object@priceOutside) +
                          (minSigma-1)*log(sharesNests)


              names(meanval)   <- object@labels


              object@slopes      <- list(alpha=minAlpha,sigma=minSigmaOut,meanval=meanval)
              object@shareInside <- 1-shareOut

              return(object)

          }

          )


logit.nests.alm <- function(prices,shares,margins,
                            ownerPre,ownerPost,
                            nests=rep(1,length(shares)),
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
       maxNests  <- nMargins - 2

       if(nNestParm > maxNests){
           stop("Additional margins must be supplied in order to calibrate nesting parameters")
       }

       if(missing(parmsStart)){

           nNests <- nlevels(nests)
           parmsStart <- runif(nNests+2) # nesting parameter values are assumed to be between 0 and 1
           parmsStart[1] <- -1* parmsStart[1] # price coefficient is assumed to be negative

           if(constraint){parmsStart <- parmsStart[1:3]}
       }


       if(constraint && length(parmsStart)!=3){
           stop("when 'constraint' is TRUE, 'parmsStart' must be a vector of length 3")
       }
       else if(!constraint && nNestParm + 2 != length(parmsStart)){
           stop("when 'constraint' is FALSE, 'parmsStart' must be a vector of length ",nNestParm + 2)

       }

       ## Create Logit  container to store relevant data
       result <- new("LogitNestsALM",prices=prices, shares=shares,
                     margins=margins,
                     ownerPre=ownerPre,
                     ownerPost=ownerPost,
                     nests=nests,
                     mcDelta=mcDelta,
                     subset=subset,
                     priceOutside=priceOutside,
                     priceStart=priceStart,
                     shareInside=sum(shares),
                     parmsStart=parmsStart,
                     constraint=constraint,
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


