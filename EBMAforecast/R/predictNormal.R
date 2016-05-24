#' @rdname EBMApredict
setMethod(f="prediction",
          signature(EBMAmodel="FDatFitNormal"),
          definition=function(EBMAmodel, 
                              Predictions,
                              Outcome=c(),
                              ...)
          {
            nDraws = 1
            if(is.matrix(Predictions)==TRUE){
              Predictions = array(Predictions,dim=c(dim(Predictions),nDraws))
            }
            #extract variables from EBMAmodel
            modelParams = EBMAmodel@modelParams
            W = EBMAmodel@modelWeights
            modelNames = EBMAmodel@modelNames
            nMods = length(W)
            nObsTest = dim(Predictions)[1]
            nMod = length(W)
            .models = EBMAmodel@modelResults
            predType = EBMAmodel@predType
            useModelParams = EBMAmodel@useModelParams
            #if(sum(EBMAmodel@modelParams[1,,])==0 & sum(EBMAmodel@modelParams[2,,])==nMods){
            #  useModelParams = FALSE
            #}
            sigma2 = EBMAmodel@variance
            .sdVec <- rep(sqrt(sigma2), nMod) 
            

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
     
            if(useModelParams==TRUE){ 
              predTestAdj <- array(NA, dim=c(nObsTest, nMod, nDraws))
                for (k in 1:nMod){
                 for (j in 1:nDraws){
                   predTestAdj[,k,j] <- .predictTest(Predictions[,k,j], i=k)
                   }
                }
              } 
            if(useModelParams==FALSE){predTestAdj <- Predictions}
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
                .ebmaMedian(.W, .x, ..sdVec)
              }
              bmaPredTest <- array(aaply(.flatPredsTest, 1, .altQBMAnormal),  dim=c(nObsTest, 1,nDraws))
              bmaPredTest[,,-1] <- NA
            }
             
            test <- abind(bmaPredTest, Predictions, along=2);  colnames(test) <- c("EBMA", modelNames)
            if(is.null(Outcome)==TRUE){Outcome = rep(numeric(0),dim(test)[1])}
                                       


            new("FDatFitNormal",
                predTest=test,
                outcomeTest= Outcome,
                modelNames=modelNames,
                modelWeights=W,
                call=match.call()
                )
          }
          )

