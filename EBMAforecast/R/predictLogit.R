#' @rdname EBMApredict
setGeneric(name="prediction",
           def=function( EBMAmodel, 
                         Predictions,
                         Outcome,
                         ...)
           {standardGeneric("prediction")}
           )

#' @importFrom plyr alply aaply laply
#' @rdname EBMApredict
setMethod(f="prediction",
          signature(EBMAmodel="FDatFitLogit"),
          definition=function(EBMAmodel, 
                              Predictions,
                              Outcome,
                              ...)
          {

          #Outcome <- matrix(Outcome)    
             nDraws <- dim(predCalibration)[3]
            if(is.matrix(Predictions)==TRUE){
              Predictions = array(Predictions,dim=c(dim(Predictions),nDraws))
          }


            
           
            
            #extract variables and observations from EBMAmodel
            predCalibration <- slot(EBMAmodel, "predCalibration")
            predCalibration <- predCalibration[,which(names(predCalibration[1,,1])!="EBMA"),1:nDraws,drop=FALSE]
            outcomeCalibration <- slot(EBMAmodel, "outcomeCalibration")
           
            W = EBMAmodel@modelWeights
            modelNames = EBMAmodel@modelNames
            nMods = length(W)
            exp = EBMAmodel@exp
            nObsTest = dim(Predictions)[1]
            useModelParams = EBMAmodel@useModelParams


              ## Set constants
            nMod <-  ncol(predCalibration)
            nObsCal <- nrow(predCalibration); nObsTest <- nrow(Predictions)
            ZERO<-1e-4
            
            #if(sum(EBMAmodel@modelParams[1,,])==0 & sum(EBMAmodel@modelParams[2,,])==nMods){
             # useModelParams = FALSE
            #}
            
              .predictCal <- function(x){
              .rawPred <- predict(x, type="response")
              .outPred <- rep(NA, nObsCal)
              .outPred[as.numeric(names(.rawPred))] <- .rawPred
              return(.outPred)
            }

            .makeAdj <- function(x){
              .adjPred <- qlogis(x)
              .negative <- .adjPred<0
              .pos <- .adjPred>1
              .adjPred <- ((1+abs(.adjPred))^(1/exp))-1
              .miss <- is.na(.adjPred)
              .negative[.miss] <- FALSE
              .adjPred[.negative] <- .adjPred[.negative]*(-1)
              #.adjPred[.pos] <- NA
              .adjPred[.miss] <- NA
              .adjPred
            }
            
            .modelFitter <- function(preds){
              .adjPred <- .makeAdj(preds)
              .thisModel <- glm(outcomeCalibration~.adjPred, family=binomial(link = "logit"))
              if (!.thisModel$converged){stop("One or more of the component logistic regressions failed to converge.  This may indicate perfect separtion or some other problem.  Try the useModelParams=FALSE option.")}
              return(.thisModel)
            }

            .predictTest <- function(x, i){
              .models[[i]]
              temp <- matrix(x,ncol=1)
              .rawPred <- predict(.models[[i]], newdata=data.frame(.adjPred=x), type="response")
              .outPred <- rep(NA, nObsTest)
              .outPred[as.numeric(names(.rawPred))] <- .rawPred
              return(.outPred)
            }

            
            ## Fit Models
            if(useModelParams){
              .models <- alply(predCalibration, 2:3, .fun=.modelFitter)
            }

             ## Extract needed info
            if(nDraws==1 & useModelParams==TRUE){
              predCalibrationAdj <- aperm(array(laply(.models, .predictCal), dim=c(nMod, nObsCal, nDraws)), c(2,1,3))
              dim(predCalibrationAdj)
              array(laply(.models, coefficients), dim=c(nMod, 2, nDraws))
              modelParams <- aperm(array(laply(.models, coefficients), dim=c(nMod, 2, nDraws)), c(2,1,3))
            }

            if(nDraws>1 & useModelParams==TRUE){ # This code is in development for exchangeability
              predCalibrationAdj <- aperm(aaply(.models, 1:2, .predictCal), c(3,1,2))
              modelParams <- aperm(aaply(.models, 1:2, coefficients), c(3,1,2))
            }
            if(useModelParams==FALSE){
              .adjPred <- .makeAdj(predCalibration)
              .adjPred[outcomeCalibration==0,,1]<-(1-plogis(.adjPred[outcomeCalibration==0,,1]))
              .adjPred[outcomeCalibration==1,,1]<-(plogis(.adjPred[outcomeCalibration==1,,1]))
              predCalibrationAdj <- .adjPred
              modelParams <- array(c(0,1), dim=c(2,nMod,nDraws))
            }

            dimnames(modelParams) <- list(c("Constant", "Predictor"), modelNames, 1:nDraws)
            dimnames(predCalibrationAdj) <- list(1:nObsCal, modelNames, 1:nDraws)

            if(useModelParams==TRUE){
                .adjPred <- .makeAdj(Predictions)
                predTestAdj <- array(NA, dim=c(nObsTest, nMod, nDraws))
                for (k in 1:nMod){
                  for (j in 1:nDraws){
                    predTestAdj[,k,j] <- .predictTest(.adjPred[,k,j], i=k)
                  }
                }
              } 
              if(useModelParams==FALSE & is.null(Outcome)==FALSE){
                .adjPred <- .makeAdj(Predictions)
                .adjPred[Outcome==0,,1]<-(1-plogis(.adjPred[Outcome==0,,1]))
              	.adjPred[Outcome==1,,1]<-(plogis(.adjPred[Outcome==1,,1]))
                predTestAdj <- .adjPred
              }
              if(useModelParams==FALSE & is.null(Outcome)==TRUE){
              predTestAdj <- Predictions
              }

              .flatPredsTest <- matrix(aaply(predTestAdj, c(1,2), function(x) {mean(x, na.rm=TRUE)}), ncol=nMod)
              bmaPredTest <-array(aaply(.flatPredsTest, 1, function(x) {sum(x* W, na.rm=TRUE)}), dim=c(nObsTest, 1,nDraws))
              bmaPredTest <-  bmaPredTest/array(t(W%*%t(1*!is.na(.flatPredsTest))), dim=c(nObsTest, 1, nDraws))
              bmaPredTest[,,-1] <- NA

              test <- abind(bmaPredTest, Predictions, along=2);  colnames(test) <- c("EBMA", modelNames)
              if(is.null(Outcome)==TRUE){Outcome = rep(numeric(0),nObsTest)}
                new("FDatFitLogit",
                predTest=test,
                outcomeTest= Outcome,
                modelNames=modelNames,
                modelWeights=W,
                modelParams=modelParams,
                call=match.call()
                )
          }
          )


