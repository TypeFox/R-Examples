#' Predict method for hybridEnsemble objects
#'
#' Prediction of new data using a hybridEnsemble model.
#' 
#' @param object An object of class hybridEnsemble created by the function  \code{hybridEnsemble}
#' @param newdata A data frame with the same predictors as in the training data
#' @param verbose TRUE or FALSE. Should information be printed to the screen
#' @param ... Not currently used
#' 
#' @examples
#' 
#' data(Credit)
#' 
#' \dontrun{
#' hE <-hybridEnsemble(x=Credit[1:100,names(Credit) != 'Response'],
#'                     y=Credit$Response[1:100],
#'                     RF.ntree=50,
#'                     AB.iter=50,
#'                     NN.size=5,
#'                     NN.decay=0,
#'                     SV.gamma = 2^-15,
#'                     SV.cost = 2^-5,
#'                     SV.degree=2,
#'                     SV.kernel='radial')
#' 
#' predictions <- predict(hE, newdata=Credit[1:100,names(Credit) != 'Response'])
#' }
#' 
#' @references Ballings, M., Vercamer, D., Van den Poel, D., Hybrid Ensemble: Many Ensembles is Better Than One, Forthcoming.
#' @seealso \code{\link{hybridEnsemble}}, \code{\link{CVhybridEnsemble}}, \code{\link{importance.hybridEnsemble}}, \code{\link{plot.CVhybridEnsemble}}, \code{\link{summary.CVhybridEnsemble}}
#' @return A list containing the following vectors:
#' \item{predMEAN}{Predictions combined by the simple mean}
#' \item{SB}{A label denoting the single best algorithm: RF=Random Forest, LR= Bagged Logistic Regression, AB= AdaBoost, SV=Bagged Support Vector Machines, NN=Bagged Neural Networks, KF=Kernel Factory}
#' \item{predSB}{Predictions by the single best}
#' \item{predAUTHORITY}{Predictions combined by authority}
#' ..and all the combination methods that are requested in the \code{\link{hybridEnsemble}} function.
#' @author Michel Ballings, Dauwe Vercamer, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @method predict hybridEnsemble
predict.hybridEnsemble <- function(object,newdata,verbose=FALSE, ...){
    
       newdata <- newdata[,!object$constants]
      
       #bagged logit
       predLR <- data.frame(matrix(nrow=nrow(newdata),ncol=length(object$LR)))
       for (i in 1:length(object$LR)) {
         predLR[,i] <- as.numeric(predict(object$LR[[i]],newx=data.matrix(newdata),type="response",s=object$LR.lambda))
       }
       predLR <- as.numeric(rowMeans(predLR))
       predLR <- .predict.calibrate(object=object$calibratorLR, newdata=predLR)
       
       #random forest
       predRF <- as.numeric(predict(object$RF, newdata,type="prob")[,2])
       predRF <- .predict.calibrate(object=object$calibratorRF, newdata=predRF)
    
       #adaboost
       predAB <- as.numeric(predict(object$AB, newdata,type="probs")[,2])
       predAB <- .predict.calibrate(object=object$calibratorAB, newdata=predAB)
    
       #kernel factory
       predKF <- as.numeric(predict(object$KF, newdata,type="probs"))
       predKF <- .predict.calibrate(object=object$calibratorKF, newdata=predKF)
      
       #neural networks
       predNN <- data.frame(matrix(nrow=nrow(newdata),ncol=length(object$NN)))
       for (i in 1:length(object$NN)) {
           newdatascaled <- data.frame(sapply(newdata, as.numeric))
           newdatascaled <- data.frame(sapply(newdatascaled, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))

         
           newdatascaled <- data.frame(t((t(newdatascaled) - ((object$minima[[i]] + object$maxima[[i]])/2))/((object$maxima[[i]]-object$minima[[i]])/2))) 
         
           predNN[,i] <- as.numeric(predict(object=object$NN[[i]],newdata=newdatascaled,type="raw"))  

       }
       predNN <- as.numeric(rowMeans(predNN))
       predNN <- .predict.calibrate(object=object$calibratorNN, newdata=predNN)
       
       #support vector machines
       predSV <- data.frame(matrix(nrow=nrow(newdata),ncol=length(object$SV)))
       for (i in 1:length(object$SV)) {
           predSV[,i] <- as.numeric(attr(predict(object$SV[[i]],newdata, probability=TRUE),"probabilities")[,2])
           
       }      
       predSV <- as.numeric(rowMeans(predSV))
       predSV <- .predict.calibrate(object=object$calibratorSV, newdata=predSV)
       
       #rotation forest
       predRoF <- as.numeric(predict(object$RoF,newdata))
       predRoF <- .predict.calibrate(object=object$calibratorRoF, newdata=predRoF)
       
       #k-nearest neighbors
  

       newdata_KNN <- data.frame(sapply(newdata, as.numeric))
       newdata_KNN <- data.frame(sapply(newdata_KNN, function(x) if(length(unique(x))==2 && min(x)==1) x-1 else x))
         
       predKNN <- data.frame(matrix(nrow=nrow(newdata_KNN),ncol=object$KNN.size))
       for (i in 1:object$KNN.size){
          ind <- sample(1:nrow(object$x_KNN),size=round(nrow(object$x_KNN)), replace=TRUE)
          #retrieve the indicators of the k nearest neighbors of the query data 
          indicatorsKNN <- as.integer(knnx.index(data=object$x_KNN, query=newdata_KNN, k=object$KNN.K))
          #retrieve the actual y from the tarining set
          predKNNoptimal <- as.integer(as.character(object$y_KNN[indicatorsKNN]))
          #if k > 1 than we take the proportion of 1s
          predKNN[,i] <- rowMeans(data.frame(matrix(data=predKNNoptimal,ncol=object$KNN.K,nrow=nrow(newdata_KNN))))
       }
        
       predKNN <- rowMeans(predKNN)
       predKNN <- .predict.calibrate(object=object$calibratorKNN, newdata=predKNN)
  
  
       #####     
  
       predictions <- data.frame(LR=predLR,RF=predRF,AB=predAB,KF=predKF,NN=predNN,SV=predSV,RoF=predRoF,KN=predKNN)
    
      
       result <- list()

       if (tolower('rbga') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Genetic Algorithm \n') 
         result$predRBGA <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsRBGA))))
       }
       if (tolower('DEopt') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Differential Evolutionary Algorithm \n')     
         result$predDEOPT <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsDEOPT))))
       }
       if (tolower('GenSA') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Generalized Simulated Annealing \n')     
         result$predGENSA <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsGENSA))))
       }
       if (tolower('malschains') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Memetic Algorithm with Local Search Chains \n')   
         result$predMALSCHAINS <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsMALSCHAINS))))
       }
       if (tolower('psoptim') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Particle Swarm Optimization \n')   
         result$predPSOPTIM <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsPSOPTIM))))
       }
       if (tolower('soma') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Self-Organising Migrating Algorithm \n') 
         result$predSOMA <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsSOMA))))
       }
       if (tolower('tabu') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Tabu Search Algorithm \n')   
         result$predTABU <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsTABU))))
       }
       
       if (tolower('LHNNLS') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Lawson-Hanson Non-negative least squares \n')    
         result$predLHNNLS <- as.numeric(crossprod(t(predictions), object$weightsLHNNLS))
       }
       
       if (tolower('GINNLS') %in% tolower(object$combine)) {
        if (verbose==TRUE) cat('   Goldfarb-Idnani Non-negative least squares \n')
        result$predGINNLS <- as.numeric(crossprod(t(predictions), object$weightsGINNLS))
       }
       
       if (tolower('NNloglik') %in% tolower(object$combine)) {
         if (verbose==TRUE) cat('   Non-negative binomial likelihood  \n')
         trimLogit <- function(x, trim=0.00001) {
           x[x < trim] <- trim
           x[x > (1-trim)] <- (1-trim)
           foo <- log(x/(1-x))
           return(foo)
         }
         
         result$predNNloglik <- as.numeric(plogis(crossprod(t(trimLogit(predictions)), object$weightsNNloglik)))
       }


       #################################################
       
       if (verbose==TRUE) cat('   Mean \n')
       result$predMEAN <- as.numeric(rowMeans(predictions))
       
       if (verbose==TRUE) cat('   Single Best \n')
       result$SB <- object$SB
       result$predSB <- predictions[,colnames(predictions) %in% object$SB]
       
       if (verbose==TRUE) cat('   Authority \n')
       result$predAUTHORITY <- as.numeric(rowSums(t(t(predictions)*as.numeric(object$weightsAUTHORITY))))
    return(result)
}




