#' @useDynLib EBMAforecast
#' @importFrom Rcpp sourceCpp


#' @importFrom plyr alply aaply
#' @rdname calibrateEnsemble
setMethod(f="fitEnsemble",
          signature(.forecastData="ForecastDataNormal"),
          definition=function(.forecastData, tol = sqrt(.Machine$double.eps),
            maxIter=1e6,
            method="EM",
            exp=numeric(),
            useModelParams = TRUE,
            predType="posteriorMedian",
            const=0, 
            W = c())
          {
            
            #check wether W is of right length and sums to 1
           if(length(W) != dim(.forecastData@predCalibration)[2] & is.null(W)==FALSE){
              stop("Vector of initial model weights must be of length of the number of predictive models included.")}  
           if(sum(W) != 1 & is.null(W)==FALSE){
              stop("Vector of initial model weights must sum to 1.")}  
            
            
          	#old EM
            # .em <- function(outcomeCalibration, prediction, RSQ, W, sigma2)
              # {

                
                # ## Step 1: Calculate the Z's
                # g<- aperm(array(aaply(.data=1:nMod,.margins=1,
                           # .fun=function(i,y, mu, sd){
                             # dnorm(y,mean=mu[,i,], sd=sd)
                           # },
                           # y=outcomeCalibration,mu=prediction, sd=sqrt(sigma2))
                          # , dim=c(nMod, nObsCal, nDraws)), c(2,1,3))
                # z.numerator<- aaply(.data=g, .margins=1, .fun=function(x){x*W})
                # z.denom <- aaply(z.numerator, 1, sum, na.rm=T)
                # Z <-aperm(array(aaply(z.numerator, 2, function(x){x/z.denom}), dim=c(nMod, nObsCal, nDraws)), c(2,1,3))
                # Z[Z < ZERO] <- 0

                # .missZ <- aaply(Z, 1, .fun=function(x) sum(!is.na(x)*1))
                # .adjConst <- const*1/.missZ
                # Z <- .adjConst + (1-const)*Z


                
                # Z[is.na(Z)] <- 0

                # ## Step 2: Calculat the W's
                # .unnormalizedW<-aaply(Z, 2, sum, na.rm = TRUE)
                # W <- .unnormalizedW
                # W <- W/sum(.unnormalizedW)
                # W[W<ZERO]<-0
                # names(W) <- modelNames
                
                # ## Step 3: Calculate sigma squared
                # sigma2<-sum(Z * RSQ, na.rm=T)/sum(Z, na.rm=T) 


                # ## Step 4: Calculate the log-likelihood
                # LL <-sum(log(z.denom))

                # out <- list(LL=LL, W=W,sigma2=sigma2)
                # return(out)
              # }
          

            .predictCal <- function(x){
              .rawPred <- predict(x)
              .outPred <- rep(NA, nObsCal)
              .outPred[as.numeric(names(.rawPred))] <- .rawPred
              return(.outPred)
            }

            
            .modelFitter <- function(preds){
              thisModel <- lm(outcomeCalibration~preds)
              return(thisModel)
            }


            .predictTest <- function(x, i){
              .rawPred <- predict(.models[[i]], newdata=data.frame(preds=x))
              .outPred <- rep(NA, nObsTest)
              .outPred[as.numeric(names(.rawPred))] <- .rawPred
              return(.outPred)
            }
            
            .ebmaMedian<-function(W, x, sdVec){
              .x <- x[!is.na(x)]
              .W <- W[!is.na(x)]
              .sdVec <- sdVec[!is.na(x)]
              
              ebmaCdf<-function(z, .x, .sdVec, .W){
                sum(.W*pnorm(z, mean=.x, sd=.sdVec))
              }
              low <- min(.x-6*.sdVec)
              up <- max(.x+6*.sdVec)
              out <- uniroot(function(z){ebmaCdf(z, .x=.x, .sdVec=.sdVec, .W=.W)-.5}
                             , lower = low
                             , upper = up
              )
              
              out$root
            }
            
            
            ##Extract data
            predCalibration <- .forecastData@predCalibration; outcomeCalibration <- .forecastData@outcomeCalibration
            predTest <- .forecastData@predTest; outcomeTest <- .forecastData@outcomeTest
            .testPeriod <- length(predTest)>0            
            modelNames <- .forecastData@modelNames
            
            ## Set constants
            nMod <-  ncol(predCalibration); nDraws <- dim(predCalibration)[3]
            nObsCal <- nrow(predCalibration); nObsTest <- nrow(predTest)
            ZERO<-1e-4
            
            ## Fit Models
            if(useModelParams==TRUE){.models <- alply(predCalibration, 2:3, .fun=.modelFitter)}

            ## Extract needed info
            if(nDraws==1 & useModelParams==TRUE){
              predCalibrationAdj <- aperm(array(laply(.models, .predictCal), dim=c(nMod, nObsCal, nDraws)), c(2,1,3))
              modelParams <- aperm(array(laply(.models, coefficients), dim=c(nMod, 2, nDraws)), c(2,1,3))
            }
            if(nDraws>1 & useModelParams==TRUE){ # This code is in development for exchangeability
              predCalibrationAdj <- aperm(aaply(.models, 1:2, .predictCal), c(3,1,2))
              modelParams <- aperm(aaply(.models, 1:2, coefficients), c(3,1,2))
            }
            if(useModelParams==FALSE){
              predCalibrationAdj <- predCalibration
              modelParams <- array(c(0,1), dim=c(2,nMod,nDraws))
            }
            calResiduals <- outcomeCalibration-predCalibrationAdj
            calResiduals2 <- calResiduals^2
            
            dimnames(modelParams) <- list(c("Constant", "Predictor"), modelNames, 1:nDraws)
            dimnames(calResiduals) <- dimnames(calResiduals2) <-dimnames(predCalibrationAdj) <- list(1:nObsCal, modelNames, 1:nDraws)

            ## Set initial values for parameters
            if(is.null(W)){
            W <- rep(1/(nMod), nMod) ; names(W) <- modelNames
            }
            sigma2<-1
            
           
            

###old EM, now run in rcpp
            # ## Run EM
            # .done <- FALSE
            # .iter <- 0
            # .emOld<-0
            # while(.done == FALSE & .iter<maxIter){
              # .thisOut <- .em(outcomeCalibration=outcomeCalibration, prediction=predCalibrationAdj, W=W, sigma2=sigma2, RSQ=calResiduals2)
              # W <- .thisOut$W
              # sigma2<-.thisOut$sigma2
              # LL <- .thisOut$LL
              # .done <- abs(.emOld-LL)/(1+abs(LL))<tol
              # .emOld <- LL
              # .iter <- .iter+1
            # }
            # if (.iter==maxIter){print("WARNING: Maximum iterations reached")}
            # W <- W*rowSums(!colSums(predCalibration, na.rm=T)==0); names(W) <- modelNames

### call to rcpp for em
			out  = emNorm(outcomeCalibration, matrix(predCalibrationAdj[,,1],ncol=nMod),matrix(calResiduals2[,,1],ncol=nMod), W, tol, maxIter, const, sigma2)
            if (out$Iterations==maxIter){print("WARNING: Maximum iterations reached")}
            W <- out$W*rowSums(!colSums(predCalibration, na.rm=T)==0); names(W) <- modelNames
            sigma2 = out$Sigma2
            LL = out$LL
            iter = out$Iterations

            ## Merge the EBMA forecasts for the calibration sample onto the predCalibration matrix
            .flatPreds <- aaply(predCalibrationAdj, c(1,2), function(x) {mean(x, na.rm=TRUE)})
            .sdVec <- rep(sqrt(sigma2), nMod) 

            if (predType=="posteriorMean"){
              bmaPred <- array(aaply(.flatPreds, 1, function(x) {sum(x* W, na.rm=TRUE)}), dim=c(nObsCal, 1,nDraws))
              bmaPred <-  bmaPred/array(t(W%*%t(1*!is.na(.flatPreds))), dim=c(nObsCal, 1, nDraws))
              bmaPred[,,-1] <- NA
            }

            if (predType=="posteriorMedian"){
              .altQBMAnormal <- function(x){
                .x <- x[!is.na(x)]
                .W <- W[!is.na(x)]
                ..sdVec <- .sdVec[!is.na(x)]
                .ebmaMedian(.W, .x, ..sdVec)
              }
             bmaPred <- array(aaply(.flatPreds, 1, .altQBMAnormal),  dim=c(nObsCal, 1,nDraws))
             bmaPred[,,-1] <- NA
            }
            cal <- abind(bmaPred, .forecastData@predCalibration, along=2); colnames(cal) <- c("EBMA", modelNames)


     
            if(.testPeriod){
              if(useModelParams==TRUE){ 
                predTestAdj <- array(NA, dim=c(nObsTest, nMod, nDraws))
                for (k in 1:nMod){
                 for (j in 1:nDraws){
                   predTestAdj[,k,j] <- .predictTest(predTest[,k,j], i=k)
                   }
                }
              } 
              if(useModelParams==FALSE){predTestAdj <- predTest}
              .flatPredsTest <- matrix(aaply(predTestAdj, c(1,2), function(x) {mean(x, na.rm=TRUE)}), ncol=nMod)
              
              if (predType=="posteriorMean"){
                bmaPredTest <-array(aaply(.flatPredsTest, 1, function(x) {sum(x* W, na.rm=TRUE)}), dim=c(nObsTest, 1,nDraws))
                bmaPredTest <-  bmaPredTest/array(t(W%*%t(1*!is.na(.flatPredsTest))), dim=c(nObsTest, 1, nDraws))
                bmaPredTest[,,-1] <- NA
              }
             
              if (predType=="posteriorMedian"){
                .altQBMAnormal <- function(x){
                  .x <- x[!is.na(x)]
                  .W <- W[!is.na(x)]
                  ..sdVec <- .sdVec[!is.na(x)]
                  .ebmaMedian( .W, .x, ..sdVec)
                }
                bmaPredTest <- array(aaply(.flatPredsTest, 1, .altQBMAnormal),  dim=c(nObsTest, 1,nDraws))
                bmaPredTest[,,-1] <- NA
              }
             
              test <- abind(bmaPredTest, .forecastData@predTest, along=2);  colnames(test) <- c("EBMA", modelNames)
            }
            if(!.testPeriod){{test <- .forecastData@predTest}}
            if(useModelParams==FALSE){.models = list()}

            new("FDatFitNormal",
                predCalibration=cal,
                outcomeCalibration=outcomeCalibration,
                predTest=test,
                outcomeTest=.forecastData@outcomeTest,
                modelNames=modelNames,
                modelWeights=W,
                useModelParams = useModelParams,
                modelParams=modelParams,
                variance=sigma2,
                logLik=LL,
                exp=exp,
                tol=tol,
                maxIter=maxIter,
                predType=predType,
                method=method,
                iter=iter,
                model="normal",
                modelResults = .models,
                call=match.call()
                )
          }
          )

