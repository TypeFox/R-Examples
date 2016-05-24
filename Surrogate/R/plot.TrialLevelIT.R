plot.TrialLevelIT <-  
  function(x, Xlab.Trial, Ylab.Trial, Main.Trial, 
           Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)),...) {
    
    Object <- x 
    if (missing(Xlab.Trial)) {Xlab.Trial <- expression(paste("Treatment effect on the Surrogate endpoint ", (alpha[i])))}
    if (missing(Ylab.Trial)) {Ylab.Trial <- expression(paste("Treatment effect on the True endpoint  ",(beta[i])))}
    if (missing(Main.Trial)) {Main.Trial <- c("Trial-level surrogacy")}
    dev.new()
    par=Par
    plot(Object$Alpha.Vector, Object$Beta.Vector, xlab=Xlab.Trial, ylab=Ylab.Trial, main=Main.Trial, ...)
    abline(lm(Object$Beta.Vector ~ Object$Alpha.Vector))
  } 
