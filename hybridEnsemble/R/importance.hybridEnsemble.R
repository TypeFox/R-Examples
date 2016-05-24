#' Importance method for hybridEnsemble objects
#'
#' Assess the importance of new data using a hybridEnsemble model. The importance is computed as follows. For each variable, compute the AUC of the model before permuting that variable and after. Next, subtract the latter from the former. This is called the decrease in AUC. If CV is greater than one, the mean is taken from all runs.
#' 
#' @param x An object of class hybridEnsemble created by the function \code{hybridEnsemble}
#' @param xdata A test data frame with the same predictors as in the training data
#' @param ydata A test factor of observed class labels (responses) with the only allowed values \{0,1\}.
#' @param method One of 'RBGA' (Genetic Algorithm), 'DEOPT' (Differential Evolution), 'GENSA' (Generalized Simulated Annealing), 'MALSCHAINS' (Memetic Algorithm), 'PSOPTIM' (Particle Swarm), 'SOMA' (Self Organizing Migrating Algorithm), 'TABU' (Tabu Search), 'LHNNLS' (Lawson-Hanson Non-negative least squares), 'GINNLS' (Goldfarb-Idnani Non-negative least squares), 'NNloglik' (Non-negative binomial likelihood), 'MEAN' (Simple Mean), 'SB' (Single Best), 'AUTHORITY' (Authority Based method)
#' @param CV An integer indicating the number of cross-validation runs
#' @param sort TRUE or FALSE. Should the predictors be sorted with the most important ones on top? 
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
#'  importance(hE,
#'           xdata=Credit[1:100,names(Credit) != 'Response'],
#'           ydata=Credit$Response[1:100])                   
#' }
#'                      
#' 
#' @references Ballings, M., Vercamer, D., Van den Poel, D., Hybrid Ensemble: Many Ensembles is Better Than One, Forthcoming.
#' @seealso \code{\link{hybridEnsemble}}, \code{\link{predict.hybridEnsemble}}, \code{\link{CVhybridEnsemble}}, \code{\link{plot.CVhybridEnsemble}}
#' @return A data frame with two colums: the variable name and the importance of the variable.
#' @author Michel Ballings, Dauwe Vercamer, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @method importance hybridEnsemble
importance.hybridEnsemble <- function(x=NULL,xdata=NULL,ydata=NULL, method="MEAN", CV=1, sort=TRUE){
  
  method <- match.arg(toupper(method),c('RBGA','DEOPT','GENSA','MALSCHAINS','PSOPTIM','SOMA','TABU','LHNNLS','GINNLS','NNloglik','MEAN','SB','AUTHORITY'))
  
  #auc before permutation
  pred <- predict(x,xdata)[[paste("pred",method,sep='')]]
  auc <- numeric()
  auc <- performance(prediction(as.numeric(pred),ydata),"auc")@y.values[[1]]
  
  store <- list()
  

  for (i in seq_len(CV)) {
      
      #store[[i]] <- data.frame(VarName=colnames(xdata),DecreaseAUC=NA)
      store[[i]] <- numeric(ncol(xdata))
      
      #auc after permutation
      for (ii in 1:ncol(xdata)){
        xdata.after <- xdata
        xdata.after[,ii] <- sample(xdata[,ii])
  
        pred <- predict(x,xdata.after)[[i]]
        auc.after <- performance(prediction(as.numeric(pred),ydata),"auc")@y.values[[1]]
        #store[[i]][ii,2] <- auc - auc.after
        store[[i]][ii] <- auc - auc.after
      }
  
  }
  
  
  store <- data.frame(VarName=colnames(xdata),DecreaseAUC=rowMeans(do.call(cbind,store)))
  
  
  if (sort==TRUE) {
    store[order(-store[,2]),]
  } else {
    store
  }
}