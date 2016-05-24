.reg1modelBasic <-
function(lm.out, dname="mydata", digits.d=NULL, show.R=FALSE) {

  nm <- all.vars(lm.out$terms)  # names of vars in the model
  n.vars <- length(nm)
  n.pred <- n.vars - 1L
  n.obs <- nrow(lm.out$model)

  tx <- character(length = 0)

# --------------
# Basic Analysis
# --------------
 
# estimates, HTs 
  sm <- summary(lm.out)

  if (show.R) {
    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- .dash2(68)
    tx[length(tx)+1] <- "> summary(model)"
    tx[length(tx)+1] <- "> confint(model)"
    tx[length(tx)+1] <- .dash2(68)
  }
  
  # output: header
  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- "Estimated Model"
    tx[length(tx)+1] <- ""
  }

  # model coefficients
  sm1 <- sm$coefficients
  sm2 <- confint(lm.out, level=0.95) 
  smc <- cbind(sm1, sm2)

  buf <- 0
  for (i in 1:n.vars) {
    lng.lbl <- nchar(rownames(smc)[i])
    if (lng.lbl > buf) buf <- lng.lbl 
   }

  max.num <- integer(length=0)
  for (icol in 1:6) {
    max.num[icol] <- 0 
    for (i in 1:n.vars) {
      ln.nm <- nchar(as.character(trunc(smc[i,icol]))) + digits.d + 1
      if (ln.nm > max.num[icol]) max.num[icol] <- ln.nm
    }
    if (max.num[icol] < 9) max.num[icol] <- 9L 
  }

  # output: row labels
  est.lbl <- .fmtc("Estimate", max.num[1]+1)
  ste.lbl <- .fmtc("  Std Err", max.num[2]+2)
  t.lbl <-  "  t-value"
  p.lbl <-  "  p-value"
  lb.lbl <- .fmtc("Lower 95%", max.num[5]+3)
  ub.lbl <- .fmtc("Upper 95%", max.num[6]+3)
  tx[length(tx)+1] <- paste(format("", width=buf), est.lbl, ste.lbl,
                           t.lbl, p.lbl, lb.lbl, ub.lbl, sep="")

  # output: values row by row
  for (i in 1:(nrow(smc))) {
    rlb <- .fmtc(rownames(smc)[i], buf)
    est <- .fmt(smc[i,1], digits.d, max.num[1])
    ste <- .fmt(smc[i,2], digits.d, max.num[2]+1)
    tvl <- .fmt(smc[i,3], 3, 8)
    pvl <- .fmt(smc[i,4], 3, 8)
    lb <- .fmt(smc[i,5], digits.d, max.num[5])
    ub <- .fmt(smc[i,6], digits.d, max.num[6])
    tx[length(tx)+1] <- paste(rlb, est, ste, tvl, pvl, " ", lb, " ", ub)
  }

  return(list(tx=tx, estimates=smc[,1], sterrs=smc[,2], tvalues=smc[,3],
     pvalues=round(smc[,4],3), cilb=smc[,5], ciub=smc[,6]))

}
