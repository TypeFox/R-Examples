plot.FixedDiscrDiscrIT <- function(x,  Weighted=TRUE, Xlab.Trial, Ylab.Trial, Main.Trial,
                                          Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)),...){
  Object <- x 
    if (missing(Xlab.Trial)) {Xlab.Trial <- expression(paste("Treatment effect on the Surrogate endpoint ", (alpha[i])))}
    if (missing(Ylab.Trial)) {Ylab.Trial <- expression(paste("Treatment effect on the True endpoint  ",(beta[i])))}
    if (missing(Main.Trial)) {Main.Trial <- c("Trial-level surrogacy")}
    dev.new()
    par=Par
    if (Weighted==TRUE){
      plot(Object$Trial.Spec.Results$Treatment.S, Object$Trial.Spec.Results$Treatment.T, cex=(Object$Trial.Spec.Results$Obs.per.trial/max(Object$Trial.Spec.Results$Obs.per.trial))*8, xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial,...)
      abline(lm(Object$Trial.Spec.Results$Treatment.T ~ Object$Trial.Spec.Results$Treatment.S))}
    
    if (Weighted==FALSE){
      plot(Object$Trial.Spec.Results$Treatment.S, Object$Trial.Spec.Results$Treatment.T, xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial, ...)
      abline(lm(Object$Trial.Spec.Results$Treatment.T ~ Object$Trial.Spec.Results$Treatment.S))}
  }