plot.TwoStageSurvSurv <- function(x, Weighted=TRUE, xlab, ylab, main,  
  Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)),...) {
 
  Object <- x 
    if (missing(xlab)) {Xlab.Trial <- expression(paste("Treatment effect on the Surrogate endpoint ", (alpha[i])))}
    if (missing(ylab)) {Ylab.Trial <- expression(paste("Treatment effect on the True endpoint  ",(beta[i])))}
    if (missing(main)) {Main.Trial <- c("Trial-level surrogacy")}
    dev.new()
    par=Par
    if (Weighted==TRUE){
      plot(Object$Results.Stage.1$LogHazard_True, Object$Results.Stage.1$LogHazard_Surr, 
           cex=(Object$Results.Stage.1$Trial.Size/max(Object$Results.Stage.1$Trial.Size))*8, 
           xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial, ...)
      abline(lm(Object$Results.Stage.1$LogHazard_True ~ Object$Results.Stage.1$LogHazard_Surr))}
    
    if (Weighted==FALSE){
      plot(Object$Results.Stage.1$LogHazard_True, Object$Results.Stage.1$LogHazard_Surr, 
           cex=1, xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial, ...)
      abline(lm(Object$Results.Stage.1$LogHazard_True ~ Object$Results.Stage.1$LogHazard_Surr))}
  
} 