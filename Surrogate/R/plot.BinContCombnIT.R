plot.FixedBinBinIT <- plot.FixedBinContIT <- plot.FixedContBinIT <- function(x, Trial.Level=TRUE, Weighted=TRUE, Indiv.Level.By.Trial=TRUE, 
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
      plot(Object$Trial.Spec.Results$Treatment.S, Object$Trial.Spec.Results$Treatment.T, 
           cex=(Object$Trial.Spec.Results$Obs.per.trial/max(Object$Trial.Spec.Results$Obs.per.trial))*8, 
           xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial,...)
      abline(lm(Object$Trial.Spec.Results$Treatment.T ~ Object$Trial.Spec.Results$Treatment.S))}
    
    if (Weighted==FALSE){
      plot(Object$Trial.Spec.Results$Treatment.S, Object$Trial.Spec.Results$Treatment.T, 
           xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial, ...)
      abline(lm(Object$Trial.Spec.Results$Treatment.T ~ Object$Trial.Spec.Results$Treatment.S))}
  }
  
  if (Indiv.Level.By.Trial==TRUE){ 
    if (missing(Xlab.Indiv)) {Xlab.Indiv <- expression(paste("R"[h.ind]^{2}))}
    if (missing(Ylab.Indiv)) {Ylab.Indiv <- "Trial"}
    if (missing(Main.Indiv)) {Main.Indiv <- c("Individual-level surrogacy")}
    dev.new()
    par=Par
    plot(y=(1:dim(Object$R2h.Ind.By.Trial)[1]), x=Object$R2h.Ind.By.Trial[,2], yaxt = "n",
         xlab=Xlab.Indiv, ylab=Ylab.Indiv, main=Main.Indiv, xlim=c(0, 1), ...)
    axis(2, at=(1:dim(Object$R2h.Ind.By.Trial)[1]), labels=Object$R2h.Ind.By.Trial[,1])
    for (i in 1: dim(Object$R2h.Ind.By.Trial)[1]){
      segments(y0 = i, x0 = Object$R2h.Ind.By.Trial[i,3], y1 = i, x1 = Object$R2h.Ind.By.Trial[i,4], lty=2)
    }
    
  }    
} 