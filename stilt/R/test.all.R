test.all <-
function(final.emul, t.plot) {
  
    # PRELIMINARIES #!+
    #==================
    p.par <- final.emul$p
    n.par <- final.emul$n


    # PERFORM CROSS-VALIDATION #!+
    #=========================
    # PRELIMINARIES #!+
    # model.out.all, emul.out.all => Matrices of model output and emulator predictions
    # [row, col] = [tested run index, time index] #!+
    model.out.all   <- matrix(NA, nrow=p.par, ncol=n.par)
    emul.out.all    <- matrix(NA, nrow=p.par, ncol=n.par)
    # Emulator standard deviation and normalized (by model range at each time) std
    # (same format) #!+
    emul.std.all    <- matrix(NA, nrow=p.par, ncol=n.par)
    emul.std.all.nm <- matrix(NA, nrow=p.par, ncol=n.par)
    # Model range for each time #!+
    model.range     <- matrix(NA, nrow=1,   ncol=n.par)


    # PERFORM CROSS-VALIDATION #!+
    for (run.out in 1:p.par) {
       csv.out <- test.csv(final.emul, 1, FALSE, NULL, run.out, FALSE) 
       model.out.all[run.out,] <- csv.out$model.out.test
       emul.out.all[run.out,]  <- csv.out$emul.out.test
       emul.std.all[run.out,]  <- csv.out$emul.std.test
    }

    # CALCULATE NORMALIZED ERRORS #!+
    # Model extremes for each time index #!+
    model.max               <- apply(model.out.all, 2, max) 
    model.min               <- apply(model.out.all, 2, min) 
    model.range             <- model.max - model.min
    stopifnot(model.range[t.plot] != 0) #Exception for zero range

    # Calculate normalized errors #!+
    model.range.mat         <- matrix(model.range, nrow=p.par, ncol=n.par, byrow=TRUE) 
    emul.std.all.nm         <- ((model.out.all - emul.out.all)/model.range.mat)*100 


    # PLOT #!+
    #======
    # Preliminaries #!+
    mfrow0 <- par("mfrow")
    on.exit(par(mfrow=mfrow0))
    par(mfrow=c(1,2))

    # x/y plot #w
    plot(model.out.all[,t.plot], emul.out.all[,t.plot], xlab="True Output",
       ylab="Predicted Output", pch=19, cex.axis=1, cex.lab=1)
    lines(c(min(1.5*model.out.all), max(1.5*model.out.all)), c(min(1.5*model.out.all),
            max(1.5*model.out.all)), col="blue") 

    # Relative errors (%) #!+
    plot(emul.std.all.nm[,t.plot], xlab="Excluded Run Number", ylab="Prediction Error (% of total range)",
         col="black", pch=19, cex.axis=1, cex.lab=1)
    lines(c(-p.par, 1.5*p.par), c(0,0), lty=2, col="black")

}
