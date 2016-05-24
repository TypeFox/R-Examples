#' Plot the performance of the cross-validated Hybrid Ensemble
#'
#' This function plots the averaged ROC curve per combination method or the median predictive performance (Area under the ROC, sensitivity or specificity curve depending on what was used in the \code{CVhybridEnsemble} function).
#' 
#' @param x An object of class CVhybridEnsemble
#' @param y Not used
#' @param ROCcurve TRUE or FALSE. Should the ROC curve be plotted or the median predictive performances?
#' @param averaging For the ROC curve: "threshold" averaging, "horizontal" averaging, or "vertical" averaging.
#' @param ... Not used
#' @details In the output: 'RBGA' (Genetic Algorithm), 'DEOPT' (Differential Evolution), 'GENSA' (Generalized Simulated Annealing), 'MALSCHAINS' (Memetic Algorithm), 'PSOPTIM' (Particle Swarm), 'SOMA' (Self Organizing Migrating Algorithm), 'TABU' (Tabue Search), 'LHNNLS' (Lawson-Hanson Non-negative least squares), 'GINNLS' (Goldfarb-Idnani Non-negative least squares), 'NNloglik' (Non-negative binomial likelihood), 'MEAN' (Simple Mean), 'SB' (Single Best), 'AUTHORITY' (Authority Based method). SB names denote the single best for all cross-validation runs: RF= Random Forest, SV= Bagged Support Vector Machines, KF= Kernel Factory, AB=AdaBoost, LR=Bagged Logistic Regression, NN=Bagged Neural Networks, RoF= Rotation Forest, KN= K-Nearest Neighbors.
#' @examples
#' 
#' 
#' data(Credit)
#' 
#' \dontrun{
#' CVhE <- CVhybridEnsemble(x=Credit[1:200,names(Credit) != 'Response'],
#'                     y=Credit$Response[1:200],
#'                     verbose=TRUE,
#'                     RF.ntree=50,
#'                     KF.rp=1,
#'                     AB.iter=50,
#'                     NN.size=5,
#'                     NN.decay=0,
#'                     SV.gamma = 2^-15,
#'                     SV.cost = 2^-5,
#'                     SV.degree=2,
#'                     SV.kernel='radial')
#' 
#' plot(x=CVhE,ROCcurve= FALSE)
#' plot(x=CVhE,ROCcurve= TRUE)
#' }
#' @references Ballings, M., Vercamer, D., Van den Poel, D., Hybrid Ensemble: Many Ensembles is Better Than One, Forthcoming.
#' @seealso \code{\link{hybridEnsemble}}, \code{\link{predict.hybridEnsemble}}, \code{\link{importance.hybridEnsemble}}, \code{\link{CVhybridEnsemble}}, \code{\link{summary.CVhybridEnsemble}}
#' @author Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
#' @method plot CVhybridEnsemble
plot.CVhybridEnsemble <- function(x,y=NULL, ROCcurve= FALSE, averaging="threshold", ...) {

if (ROCcurve==TRUE){
  
  colors <- c("blue","red","green","violet","purple","cornflowerblue","orange","brown","black")[1:length(x)]
  
  for (i in 1:(length(x)-1)) {
      
    PREDS <- list()
    for (ii in 1:10) PREDS[[ii]] <- x[[i]]$predictions[[ii]]$predicted
  
    LABS <- list()
    for (ii in 1:10) LABS[[ii]] <- x[[i]]$predictions[[ii]]$response
  
    if (i==1) {
      plot(performance(prediction(PREDS, LABS), "tpr", "fpr"), col=colors[i], avg=averaging)
    }else{
      plot(performance(prediction(PREDS, LABS), "tpr", "fpr"), col=colors[i], avg=averaging, add=TRUE)
    }
    
    lines(x=seq(0,1,0.1),y=seq(0,1,0.1),lty=3)
    
    xlabels <- names(x)[1:(length(x)-1)]
    lab <- character()

    for (i in 1:length(unlist(x$SB$SBname)))  {
   
    if (i<length(unlist(x$SB$SBname))) {
      lab[i] <- paste(substr(unlist(x$SB$SBname)[i],1,2),",",sep="") 
    } else {
      lab[i] <- paste(substr(unlist(x$SB$SBname)[i],1,2)) 
    }
  }
    
    xlabels[length(xlabels)] <- paste(xlabels[length(xlabels)],":",paste(lab,collapse=""),sep="")
    
    legend("bottomright", legend=xlabels,col=colors,lty=1,lwd=3,cex=1)
  }
      
}else {
  dat <- numeric()
  for (i in 1:(length(x)-1)) dat[i] <- round(x[[i]]$median,4)
  xlabels <- names(x)[1:(length(x)-1)]

  lab <- character()

      for (i in 1:length(unlist(x$SB$SBname)))  {
   
    if (i<length(unlist(x$SB$SBname))) {
      lab[i] <- paste(substr(unlist(x$SB$SBname)[i],1,2),",",sep="") 
    } else {
      lab[i] <- paste(substr(unlist(x$SB$SBname)[i],1,2)) 
    }
  }
  xlabels[length(xlabels)] <- paste(xlabels[length(xlabels)],"\n",paste(lab,collapse=""))
  
  v <-  barplot(dat,
        ylim=c(min(dat)-0.02,max(dat)+0.02),
        names.arg=xlabels,
        cex.names=0.8,xpd=FALSE, main="5x2f CV", ylab = "Median predictive performance")


  text(x=v, y=dat+0.005, labels=dat)
}


}