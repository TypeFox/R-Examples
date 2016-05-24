plot.Single.Trial.RE.AA <- function(x, Trial.Level=TRUE, Indiv.Level=TRUE, Xlab.Indiv, Ylab.Indiv, Xlab.Trial, Ylab.Trial, 
                                    Main.Trial, Main.Indiv, Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), ...){
  
  Object <- x  
 
  if (Trial.Level==TRUE){ 
    if (missing(Xlab.Trial)) {Xlab.Trial <- expression(paste("Treatment effect on the surrogate endpoint ", (alpha)))}
    if (missing(Ylab.Trial)) {Ylab.Trial <- expression(paste("Treatment effect on the true endpoint  ",(beta)))}
    if (missing(Main.Trial)) {Main.Trial <- expression(paste("Relative Effect (RE)"))}
    dev.new()
    par=Par
    plot(Object$Alpha$Alpha, Object$Beta$Beta, xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial, ...)
    abline(a=0, b=Object$RE.Delta$RE)
  }
  
  if (Indiv.Level==TRUE){ 
    if (missing(Xlab.Indiv)) {Xlab.Indiv <- expression(paste("Residuals for the surrogate endpoint ", (epsilon[Sj])))}
    if (missing(Ylab.Indiv)) {Ylab.Indiv <- expression(paste("Residuals for the true endpoint  ", (epsilon[Tj])))}
    if (missing(Main.Indiv)) {Main.Indiv <- expression(paste("Adjusted Assocation ", (rho[Z])))}
    dev.new()
    par=Par
    plot(x=Object$Residuals$Surr, y=Object$Residuals$True, xlab=Xlab.Indiv, ylab=Ylab.Indiv, main=Main.Indiv, ...)
    abline(lm(Object$Residuals$True ~ -1 + Object$Residuals$Surr))
  }    
} 
