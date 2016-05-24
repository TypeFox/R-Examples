plot.MixedContContIT <- plot.FixedContContIT <- function(x, Trial.Level=TRUE, Weighted=TRUE, Indiv.Level=TRUE, 
                                                           Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, Main.Trial, Main.Indiv, 
                                                         Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)),...){
  Object <- x 
  if (Trial.Level==TRUE){ 
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
  
  if (Indiv.Level==TRUE){ 
    if (missing(Xlab.Indiv)) {Xlab.Indiv <- expression(paste("Residuals for the Surrogate endpoint ", (epsilon[Sij])))}
    if (missing(Ylab.Indiv)) {Ylab.Indiv <- expression(paste("Residuals for the True endpoint  ", (epsilon[Tij])))}
    if (missing(Main.Indiv)) {Main.Indiv <- c("Individual-level surrogacy")}
    dev.new()
    par=Par
    plot(x=Object$Residuals$Residuals.Model.S, y=Object$Residuals$Residuals.Model.T, xlab=Xlab.Indiv, ylab=Ylab.Indiv, main=Main.Indiv, ...)
    abline(lm(Object$Residuals$Residuals.Model.T ~ -1 + Object$Residuals$Residuals.Model.S))
  }    
} 

