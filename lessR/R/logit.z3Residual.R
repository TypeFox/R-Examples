.logit3Residual <-
function(lm.out, nm, mydata,
         n.vars, n.pred, n.obs, n.keep, digits.d, pre, line,
         res.sort, res.rows, cooks.cut) {
  
    cat( "\n\n\n", "  ANALYSIS OF RESIDUALS AND INFLUENCE", "\n")

    cat("Data, Fitted, Residual, Studentized Residual, Dffits, Cook's Distance\n")
    if (res.sort == "cooks") cat("   [sorted by Cook's Distance]\n")
    if (res.sort == "rstudent")  
      cat("   [sorted by Studentized Residual, ignoring + or - sign]\n")
   if (res.sort == "dffits")  
      cat("   [sorted by dffits, ignoring + or - sign]\n")
    txt <- "cases (rows) of data]"
    cat("   [res.rows = ", res.rows, " out of ", n.keep, " ", txt, sep="", "\n")
    .dash(68)

    fit <- fitted(lm.out)
    res <- residuals(lm.out, type="response")
    cook <- cooks.distance(lm.out)
    
    out <- cbind(fit, res, rstudent(lm.out), dffits(lm.out), cook)
    out <- cbind(lm.out$model[c(nm[seq(2,n.vars)],nm[1])],out)
    out <- data.frame(out)
    names(out)[n.vars+1] <- "fitted"
    names(out)[n.vars+2] <- "residual"
    names(out)[n.vars+3] <- "rstudent"
    names(out)[n.vars+4] <- "dffits"
    names(out)[n.vars+5] <- "cooks"
    if (res.sort != "off") {
      if (res.sort == "cooks") o <- order(out$cooks, decreasing=TRUE)
      if (res.sort == "rstudent") o <- order(abs(out$rstudent),
        decreasing=TRUE)
      if (res.sort == "dffits") o <- order(abs(out$dffits),
        decreasing=TRUE)
      out <- out[o,]
    }
    print(out[1:res.rows,], digits=digits.d)
    rm(out)

}
