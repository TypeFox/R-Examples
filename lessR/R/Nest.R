Nest <-
function(y, nested.model, full.model, method=c("ls", "logit"),
         data=mydata, digits.d=NULL) {

  method <- match.arg(method)

  my.vars <- as.list(seq_along(data))
  names(my.vars) <- names(data)


  # get response variable
  rv <- eval(substitute(y), envir=my.vars, enclos=parent.frame())
  y.name <-  names(my.vars)[rv]

  if (is.null(digits.d)) digits.d <- .getdigits(data[,y.name], 2)
  options(digits.d=digits.d) 

  if (method == "logit") {
    is.bin <- TRUE
    if (is.factor(data[,y.name])) { 
       if (nlevels(data[,y.name]) != 2) is.bin  <- FALSE
    }
    else {
      for (i in 1:nrow(data))
        if (!is.na(data[i,y.name]))
          if (data[i,y.name]!=0 && data[i,y.name]!=1) is.bin <- FALSE
     }
    if (!is.bin) { 
      cat("\n"); stop(call.=FALSE, "\n","------\n",
        "Response variable: ", y.name, "\n",
        "If numeric, can only have values of 0 or 1.\n",
        "If a factor, can only have two levels.\n\n")
    }
  }

  # get nested model
  nest.vars <- eval(substitute(nested.model), envir=my.vars, enclos=parent.frame())
  n.pred <- length(nest.vars)
  x.name <- character(length=n.pred)
  for (ivar in 1:n.pred) {
    x.name[ivar] <- names(my.vars)[nest.vars[ivar]]
  }
  xn.preds <- x.name[1]
  if (n.pred > 1)
    for(ivar in 2:n.pred) xn.preds <- paste(xn.preds, "+", x.name[ivar])

  n.formula <- as.formula(paste(y.name, "~", xn.preds))


  # get full model and analyze
  full.vars <- eval(substitute(full.model), envir=my.vars, enclos=parent.frame())
  full.vars <- union(nest.vars, full.vars)
  n.pred <- length(full.vars)

  x.name <- character(length=n.pred)

  for (ivar in 1:n.pred) x.name[ivar] <- names(my.vars)[full.vars[ivar]]

  xf.preds <- x.name[1]
  if (n.pred > 1)
    for(ivar in 2:n.pred) xf.preds <- paste(xf.preds, "+", x.name[ivar])

  f.formula <- as.formula(paste(y.name, "~", xf.preds))


  # do the model comparison

  if (method=="ls") {
    lm.full <- lm(f.formula, data=data)
    lm.nest <- lm(n.formula, data=lm.full$model)
    av <- anova(lm.nest, lm.full)

    tx <- character(length = 0)

    models <- attr(av, "heading")[2]
    models <- sub("Model 1", "Reduced Model", models)
    models <- sub("Model 2", "Full Model   ", models)
    tx[length(tx)+1] <- models
    txtmdl <- tx

    tx <- character(length = 0)

    if (is.null(options()$knitr.in.progress)) {
      tx[length(tx)+1] <- "Analysis of Variance"
      tx[length(tx)+1] <- ""
    }

    max.c1 <- nchar("Residual")

    mr.df <- av$Df[2]; mr.ss <- av$`Sum of Sq`[2]; mr.ms <- mr.ss/mr.df
    mr.f <- av$F[2]; mr.p <- av$`Pr(>F)`[2]

    mf.df <- av$Res.Df[2]; mf.ss <- av$RSS[2]; mf.ms <- mf.ss/mf.df

    # width of columns
    d <- digits.d
    max.ln <- integer(length=0)
    max.ln[1] <- max(nchar(.fmti(mr.df)),nchar(.fmti(mf.df)),nchar("df"))
    max.ln[2] <- max(nchar(.fmt(mr.ss,d)),nchar(.fmt(mf.ss,d)),nchar("Sum Sq")) 
    max.ln[3] <- max(nchar(.fmt(mr.f,d)),nchar("Mean Sq"))
    max.ln[4] <- nchar(.fmt(mr.p,d))
    max.ln <- max.ln + 2

    df.lbl <- .fmtc("     df", max.ln[1]+1)
    SS.lbl <- .fmtc(" Sum Sq", max.ln[2]+1)
    MS.lbl <- .fmtc("Mean Sq", max.ln[3]+1)
    fv.lbl <- .fmtc("F-value", max.ln[4]+3)
    tx[length(tx)+1] <- paste(eval(format("", width=max.c1+max.ln[1]-6)),
           df.lbl, SS.lbl, MS.lbl, fv.lbl, "   p-value", sep="")

    if (n.pred > 0) {
      rlb <- .fmtc("Tested", max.c1, j="left")

      md <- .fmti(mr.df, max.ln[1]) 
      ms <- .fmt(mr.ss, digits.d, max.ln[2])
      mm <- .fmt(mr.ms, digits.d, max.ln[3])
      mf <- .fmt(mr.f, digits.d, max.ln[4]+2)
      mp <- .fmt(mr.p, 3, 9) 

      tx[length(tx)+1] <- paste(rlb, md, ms, mm, mf, mp)

      tst <- c(mr.df, mr.ss, mr.ms, mr.f, mr.p)
      names(tst) <- c("df", "ss", "ms", "fvalue", "pvalue")

    # Residuals
    rlb <- .fmtc("Residual", max.c1, j="left")
    df <- .fmti(mf.df, max.ln[1])
    SS <- .fmt(mf.ss, digits.d, max.ln[2])
    MS <- .fmt(mf.ms, digits.d, max.ln[3])

    tx[length(tx)+1] <- paste(rlb, df, SS, MS) 

    rsd <- c(mf.df, mf.ss, mf.ms)
    names(rsd) <- c("df", "ss", "ms")

    tot <- NA

    txtbl <- tx

    }
  }  # end method="ls"

  else {
    lm.full <- suppressWarnings(glm(f.formula, data=data, family="binomial"))
    lm.nest <- suppressWarnings(glm(n.formula, data=lm.full$model, family="binomial"))
    av <- anova(lm.nest, lm.full, test="Chisq")

    tx <- character(length = 0)

    models <- attr(av, "heading")[2]
    models <- sub("Model 1", "Reduced Model", models)
    models <- sub("Model 2", "Full Model   ", models)
    tx[length(tx)+1] <- models
    txtmdl <- tx

    tx <- character(length = 0)

    if (is.null(options()$knitr.in.progress)) {
      tx[length(tx)+1] <- "Analysis of Deviance"
      tx[length(tx)+1] <- ""
    }

    max.c1 <- nchar("Residual")

    mr.df <- av$Df[2]; mr.dv <- av$`Deviance`[2]; mr.p <- av$`Pr(>Chi)`[2]
    mf.df <- av$`Resid. Df`[2]; mf.dv <- av$`Resid. Dev`[2];
    mt.df <- av$`Resid. Df`[1]; mt.dv <- av$`Resid. Dev`[1];

    # width of columns
    d <- digits.d
    max.ln <- integer(length=0)
    max.ln[1] <- max(nchar(.fmti(mr.df)),nchar(.fmti(mf.df)),nchar("df"))
    max.ln[2] <- nchar(.fmt(mr.dv,d)) + 2
    max.ln <- max.ln + 2

    df.lbl <- .fmtc("     df", max.ln[1]+1)
    dv.lbl <- .fmtc(" Deviance", max.ln[2]+1)
    tx[length(tx)+1] <- paste(eval(format("", width=max.c1-max.ln[1]+4)),
           df.lbl, dv.lbl, "   p-value", sep="")

    if (n.pred > 0) {
      rlb <- .fmtc("Tested", max.c1, j="left")

      md <- .fmti(mr.df, max.ln[1]) 
      dv <- .fmt(mr.dv, digits.d, max.ln[2])
      mp <- .fmt(mr.p, 3, 9) 

      tx[length(tx)+1] <- paste(rlb, md, dv, mp)

      tst <- c(mr.df, mr.dv, mr.p)
      names(tst) <- c("df", "dev", "pvalue")

    # Residuals
    rlb <- .fmtc("Residual", max.c1, j="left")
    df <- .fmti(mf.df, max.ln[1])
    dev <- .fmt(mf.dv, digits.d, max.ln[2])

    tx[length(tx)+1] <- paste(rlb, df, dev) 

    rsd <- c(mf.df, mf.dv)
    names(rsd) <- c("df", "dev")

    # Total
    rlb <- .fmtc("Total", max.c1, j="left")
    df <- .fmti(mt.df, max.ln[1])
    dev <- .fmt(mt.dv, digits.d, max.ln[2])

    tx[length(tx)+1] <- ""
    tx[length(tx)+1] <- paste(rlb, df, dev) 

    tot <- c(mf.df, mf.dv)
    names(tot) <- c("df", "dev")

    txtbl <- tx

    }

  }  # end method="logit"


  class(txtmdl) <- "out_piece"
  class(txtbl) <- "out_piece"
  
  output <- list(
    fun.call=match.call(), out_models=txtmdl, out_anova=txtbl,
    anova_tested=tst, anova_residual=rsd, anova_total=tot
  )

  class(output) <- "out_all"

  return(output)

}
