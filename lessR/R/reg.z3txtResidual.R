.reg3txtResidual <-
function(lm.out, cook, digits.d=NULL, res.sort="cooks", res.rows=NULL, show.R=FALSE) {

  nm <- all.vars(lm.out$terms)  # names of vars in the model
  n.vars <- length(nm)
  n.pred <- n.vars - 1L
  n.keep <- nrow(lm.out$model)
  
  tx <- character(length = 0)

# ----------------------------------------------
# text output

  if (show.R) {
    tx[length(tx)+1] <- .dash2(68)
    tx[length(tx)+1] <- paste("> ", "fitted(model)", sep="", "\n")
    tx[length(tx)+1] <- paste("> ", "resid(model)", sep="", "\n")
    tx[length(tx)+1] <- paste("> ", "rstudent(model)", sep="", "\n")
    tx[length(tx)+1] <- paste("> ", "dffits(model)", sep="", "\n")
    tx[length(tx)+1] <- paste("> ", "cooks.distance(model)", sep="", "\n")
    tx[length(tx)+1] <- .dash2(68)
  }

  tx[length(tx)+1] <- "Data, Fitted, Residual, Studentized Residual, Dffits, Cook's Distance"
  if (res.sort == "cooks")
    tx[length(tx)+1] <- "   [sorted by Cook's Distance]"
  if (res.sort == "rstudent")  
    tx[length(tx)+1] <- "   [sorted by Studentized Residual, ignoring + or - sign]"
  if (res.sort == "dffits")  
    tx[length(tx)+1] <- "   [sorted by dffits, ignoring + or - sign]"
  if (res.rows < n.keep)
    txt <- "rows of data, or do res.rows=\"all\"]"
  else
    txt="]"
  tx[length(tx)+1] <- paste("   [res.rows = ", res.rows, ", out of ", n.keep, " ", txt, sep="")


  fit <- lm.out$fitted.values
  res <- lm.out$residuals
  #cook <- cooks.distance(lm.out)
  
  # text output
  out <- data.frame(fit, res, rstudent(lm.out), dffits(lm.out), cook)
  out <- data.frame(lm.out$model[nm[1]], out)
  if (n.pred > 0) out <- data.frame(lm.out$model[c(nm[seq(2,n.vars)])], out)

  #out <- data.frame(out)
  names(out)[n.vars+1] <- "fitted"
  names(out)[n.vars+2] <- "resid"
  names(out)[n.vars+3] <- "rstdnt"
  names(out)[n.vars+4] <- "dffits"
  names(out)[n.vars+5] <- "cooks"
  if (res.sort != "off") {
    if (res.sort == "cooks") {
      o <- order(out$cooks, decreasing=TRUE)
      clmn <- 0L
    }
    if (res.sort == "rstudent") {
      o <- order(abs(out$rstdnt), decreasing=TRUE)
      clmn <- 2L
    }
    if (res.sort == "dffits") {
      o <- order(abs(out$dffits), decreasing=TRUE)
      clmn <- 1L
    }
    out <- out[o,]
  }

  tx2 <- .prntbl(out[1:res.rows,], digits.d)
  for (i in 1:length(tx2)) tx[length(tx)+1] <- tx2[i]

  if (res.rows > 5  &&  res.sort != "off") {
    label.top <- numeric(length=0)
    out.top <- numeric(length=0)
    for (i in 1:5) {
      label.top[i] <- rownames(out)[i]
      out.top[i] <- out[i,(ncol(out)-clmn)]
    }
      names(out.top) <- label.top
  }
  else
    out.top <- NA

  return(list(tx=tx, resid.max=out.top))

}
