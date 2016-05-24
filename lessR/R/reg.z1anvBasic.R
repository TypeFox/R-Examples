.reg1anvBasic <-
function(lm.out, dname="mydata", digits.d=NULL, show.R=FALSE) {

  nm <- all.vars(lm.out$terms)  # names of vars in the model
  n.vars <- length(nm)
  n.pred <- n.vars - 1L
  n.obs <- nrow(lm.out$model)
  d <- digits.d

  tx <- character(length = 0)

  # ANOVA 
  smc <- anova(lm.out)
  SSE <- smc[n.vars,2]

  if (show.R) {
    tx[length(tx)+1] <- ""
    .dash2(68)
    tx[length(tx)+1] <- paste("> ","anova(model)", "\n",sep="")
    tx[length(tx)+1] <-.dash2(68)
  }

  if (is.null(options()$knitr.in.progress)) {
    tx[length(tx)+1] <- "Analysis of Variance"
    tx[length(tx)+1] <- ""
  }

  # width of column 1
  max.c1 <- max(nchar("Model"), nchar(nm[1]))
  for (i in 1:n.vars) {
    c1 <- nchar(rownames(smc)[i])
    if (c1 > max.c1) max.c1 <- c1 
   }

  # width of data columns
  max.ln <- integer(length=0)
  for (i in 1:4) {
    ln.nm <- nchar(colnames(smc)[i])
    max.ln[i] <- ln.nm + 1
    for (j in 1:nrow(smc)) {
      xjc <- .fmt(smc[j,i], d=digits.d)
      if (nchar(xjc) > max.ln[i]) max.ln[i] <- nchar(xjc)
    }
    max.ln[i] <- max.ln[i] + 1L
    if (max.ln[i] < 9L) max.ln[i] <- 9L
  }

  df.lbl <- .fmtc("     df", max.ln[1]+1)
  SS.lbl <- .fmtc(" Sum Sq", max.ln[2]+1)
  MS.lbl <- .fmtc("Mean Sq", max.ln[3]+1)
  fv.lbl <- .fmtc("F-value", max.ln[4]+1)
  tx[length(tx)+1] <- paste(eval(format("", width=max.c1-5)), df.lbl, SS.lbl,
                             MS.lbl, fv.lbl, "   p-value", sep="")

  # predictors
  if (n.pred > 1) {
    for (i in 1:(n.pred)) {
      rlb <- .fmtc(rownames(smc)[i], max.c1)
      df <- .fmti(smc[i,1], max.ln[1]-5)
      SS <- .fmt(smc[i,2], digits.d, max.ln[2])
      MS <- .fmt(smc[i,3], digits.d, max.ln[3])
      fv <- .fmt(smc[i,4], digits.d, max.ln[4])
      pv <- .fmt(smc[i,5], 3, 9)
      tx[length(tx)+1] <- paste(rlb, df, SS, MS, fv, pv)
    } 
  }

  # Model term in ANOVA table
  mdl <- NA
  if (n.pred > 0) {
    rlb <- .fmtc("Model", max.c1, j="left")

    mod.df <- 0
    for (i in 1:n.pred) mod.df <- mod.df + smc[i,1]
    md <- .fmti(mod.df, max.ln[1]-5) 

    mod.ss <- 0 
    for (i in 1:n.pred) mod.ss <- mod.ss + smc[i,2]
    ms <- .fmt(mod.ss, digits.d, max.ln[2])

    mod.ms <- mod.ss/mod.df
    mm <- .fmt(mod.ms, digits.d, max.ln[3])

    mod.f <- mod.ms/smc[n.vars, 3]
    mf <- .fmt(mod.f, digits.d, max.ln[4])

    mod.p <- pf(mod.f, mod.df, smc[n.vars,1], lower.tail=FALSE)
    mp <- .fmt(mod.p, 3, 9) 

    tx[length(tx)+1] <- paste(rlb, md, ms, mm, mf, mp)
    if (n.pred > 1) tx[length(tx)+1] <- ""

    mdl <- c(mod.df, mod.ss, mod.ms, mod.f, mod.p)
    names(mdl) <- c("df", "ss", "ms", "fvalue", "pvalue")
  }

  # Residuals
  rlb <- .fmtc(rownames(smc)[n.vars], max.c1, j="left")
  df <- .fmti(smc[n.vars,1], max.ln[1]-5)
  SS <- .fmt(smc[n.vars,2], digits.d, max.ln[2])
  MS <- .fmt(smc[n.vars,3], digits.d, max.ln[3])
  MSW <- smc[n.vars,3]
  tx[length(tx)+1] <- paste(rlb, df, SS, MS) 
  if (n.pred > 1) tx[length(tx)+1] <- ""

  rsd <- c(smc[n.vars,1], smc[n.vars,2], smc[n.vars,3])
  names(rsd) <- c("df", "ss", "ms")

  # Total
  tot <- NA
  if (n.pred > 0) {
    rlb <- .fmtc(nm[1], max.c1, j="left")

    tot.df <- mod.df + smc[n.vars,1]
    td <- .fmti(tot.df, max.ln[1]-5) 

    tot.ss <- mod.ss + smc[n.vars,2]
    ts <- .fmt(tot.ss, digits.d, max.ln[2])

    tot.ms <- tot.ss/tot.df
    tm <- .fmt(tot.ms, digits.d, max.ln[3])

    tx[length(tx)+1] <- paste(rlb, td, ts, tm) 

    tot <- c(tot.df, tot.ss, tot.ms)
    names(tot) <- c("df", "ss", "ms")
  }


  return(list(tx=tx, mdl=mdl, rsd=rsd, tot=tot, MSW=MSW))

}
