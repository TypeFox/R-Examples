search.bayesx.tex <- function(x)
{
  rval <- list()
  fam <- any(grepl("Family:", x))
  if(fam) {
    fam <- grep("Family:", x, value = TRUE)
    fam <- splitme(fam)
    nfam<- length(fam)
    j <- nfam + 1L
    ifam <- NULL
    for(i in 1L:nfam) {
      if(fam[i] == ">")
        j <- i
      if(fam[i] == "\\")
        j <- nfam + 1L
      if(i > j)
        ifam <- c(ifam, fam[i])
    }
    rval$family <- resplit(ifam)
  }
  m <- any(grepl("BAYESREG", x))
  if(m)
    m <- "MCMC"
  else {
    m <- any(grepl("remlreg", x))
    if(m)
      m <- "REML"
    else {
      m <- any(grepl("STEPWISEREG", x))
      if(m)
        m <- "STEP"
      else
        m <- NULL
    }
  }
  rval$method <- m
  obs <- any(grepl("Number of observations", x))
  if(obs) {
    obs <- grep("Number of observations:", x, value = TRUE)
    obs <- splitme(obs)
    nobs <- length(obs)
    j <- nobs + 1L
    N <- NULL
    for(i in 1L:nobs) {
      if(obs[i] == "=" || obs[i] == ">")
        j <- i
      if(obs[i]=="\\")
        j <- nobs + 1L
      if(i > j)
        N <- c(N, obs[i])
    }
    rval$N <- as.integer(resplit(N))
    if(is.na(rval$N))
      rval$N <- NULL
  }
  iter <- any(grepl("Number of Iterations:", x))
  if(iter) {
    iter <- splitme(grep("Number of Iterations:", x, value = TRUE))
    ni <- length(iter)
		j <- ni + 1L
    igrep <- NULL
    for(i in 1L:ni) {
      if(iter[i] == ">")
        j <- i
      if(iter[i] == "\\")
        j <- ni + 1L
      if(i > j)
        igrep <- c(igrep, iter[i])
    }
    rval$iterations <- as.integer(resplit(igrep))
  }
  burn <- any(grepl("Burn in:", x))
  if(burn) {
    burn <- splitme(grep("Burn in:", x, value = TRUE))
    ni <- length(burn)
    j <- ni + 1L
    burngrep <- NULL
    for(i in 1L:ni) {
      if(burn[i] == ">")
        j <- i
      if(burn[i] == "\\")
        j <- ni + 1L
      if(i > j)
        burngrep <- c(burngrep,burn[i])
    }
    rval$burnin <- as.integer(resplit(burngrep))
  }
  thin <- any(grepl("Thinning Parameter:", x))
  if(thin) {
    thin <- splitme(grep("Thinning Parameter:", x, value = TRUE))
    ni <- length(thin)
    j <- ni + 1L
    thingrep <- NULL
   for(i in 1L:ni) {
     if(thin[i] == ">")
       j <- i
     if(i > j)
       thingrep <- c(thingrep, thin[i])
    }
    rval$step <- as.integer(resplit(thingrep))
  }
  dic <- any(grepl("DIC \\\\>", x))
  if(dic) {
    dic <- splitme(grep("DIC \\\\>", x, value = TRUE))
    ni <- length(dic)
    j <- ni + 1L
    dicgrep <- NULL
    for(i in 1L:ni) {
      if(dic[i] == ">")
        j <- i
      if(dic[i] == "\\")
        j <- ni + 1L
      if(i > j)
        dicgrep <- c(dicgrep, dic[i])
    }
    rval$DIC <- as.numeric(resplit(dicgrep))
  }
  pd <- any(grepl("pD \\\\>", x))
  if(pd) {
    pd <- splitme(grep("pD \\\\>", x, value = TRUE))
    ni <- length(pd)
    j <- ni + 1L
    pdgrep <- NULL
    for(i in 1L:ni) {
      if(pd[i] == ">")
        j <- i
      if(pd[i]=="\\")
        j <- ni + 1L
      if(i > j)
        pdgrep <- c(pdgrep, pd[i])
    }
    rval$pd <- as.numeric(resplit(pdgrep))
  }
  df <- any(grep("Degrees of freedom:", x))
  if(df) {
    df <- grep("Degrees of freedom:", x, value = TRUE)
    if(length(df) > 1L)
      df <- df[2L]
    df <- splitme(df)
    ok <- FALSE
    DF <- NULL
    for(i in 1L:length(df)) {
      if(ok && df[i] != "\\")
        DF <- c(DF,df[i])
      if(df[i] == ">")
        ok <- TRUE
    }
    rval$df <- as.numeric(resplit(DF))
  }
  aic <- any(grep("AIC:", x))
  if(aic) {
    aic <- grep("AIC:", x, value = TRUE)
    aic <- splitme(aic)
    ok <- FALSE
    AIC <- NULL
    for(i in 1L:length(aic)) {
      if(ok && aic[i] != "\\")
        AIC <- c(AIC, aic[i])
      if(aic[i] == ">")
        ok <- TRUE
    }
    rval$AIC <- as.numeric(resplit(AIC))
  }
  bic <- any(grep("BIC:", x))
  if(bic) {
    bic <- grep("BIC:", x, value = TRUE)
    bic <- splitme(bic)
    ok <- FALSE
    BIC <- NULL
    for(i in 1L:length(bic)) {
      if(ok && bic[i] != "\\")
        BIC <- c(BIC, bic[i])
      if(bic[i] == ">")
        ok <- TRUE
    }
    rval$BIC <- as.numeric(resplit(BIC))
  }
  gcv <- any(grep("GCV:", x))
  if(gcv) {
    gcv <- grep("GCV:", x, value = TRUE)
    gcv <- splitme(gcv)
    ok <- FALSE
    GCV <- NULL
    for(i in 1L:length(gcv)) {
      if(ok && gcv[i] != "\\")
        GCV <- c(GCV, gcv[i])
      if(gcv[i] == ">")
        ok <- TRUE
    }
    rval$GCV <- as.numeric(resplit(GCV))
  }
  final <- any(grepl("Final Predictor:", x))
  if(final) {
    final <- i <- grep("Final Predictor:", x)
    stepfiles <- NULL
    run <- TRUE
    while(run) {
      i <- i + 1L
      if(x[i] == "\\newpage ")
        run <- FALSE
      else
        stepfiles <- c(stepfiles, x[i])
      if(i == length(x))
        run <- FALSE
    }
    final <- splitme(grep("eta", stepfiles, value = TRUE))
    grepfinal <- NULL
    for(i in 1L:length(final)) {
      check <- final[i] != "$" && final[i] != "&" && final[i] != "\\"
      if(check) {
        ok <- final[i]
        if(ok == " ")
          ok <- "!"
        grepfinal <- c(grepfinal, ok)
      }
    }
    grepfinal <- resplit(grepfinal)
    grepfinal <- sub("cdot", "*", grepfinal)
    grepfinal <- sub("eta..", "eta!", grepfinal)
    grepfinal <- sub("..gamma_0", "!gamma0", grepfinal)
    grepfinal <- splitme(grepfinal)
    stepgrep <- NULL
    ok <- TRUE
    for(i in 1L:length(grepfinal)) {
      if(grepfinal[i] == "{")
        ok <- FALSE
      if(ok) {
        take <- grepfinal[i]
        if(take == "!")
          take <- " "
        if(take == "_")
          take <- NULL
        stepgrep <- c(stepgrep, take)
      }
      if(grepfinal[i] == "}")
        ok <- TRUE	
    }
    crit <- splitme(x[grep("Final Predictor:", x) + 8L])
    fc <- fcn <- NULL
    okn <- TRUE
    okc <- FALSE
    for(i in 1L:length(crit)) {
      if(okc && crit[i] != "\\")
        fc <- c(fc, crit[i])
      if(crit[i] == "=") {
        okn <- FALSE
        okc <- TRUE
      }
      if(okn && crit[i] != "\\")
        fcn <- c(fcn, crit[i])
    }
    fc <- as.numeric(resplit(fc))
    fcn <- resplit(fcn)
    fcn <- sub(" ", "", fcn)
    eval(parse(text=paste("rval$", fcn, "<-", fc, sep = "")))
    step.final.model <- gsub("f\\(", "sx\\(", resplit(stepgrep))
    step.final.model <- gsub("gamma0", "\\(Intercept\\)", step.final.model)
    rval$step.final.model <- step.final.model
  }
  final.prop <- any(grepl("Final Properties:", x))
  if(final.prop) {
    final <- i <- grep("Final Properties:", x)
    stepfiles <- NULL
    run <- TRUE
    while(run) {
      i <- i + 1L
      if(x[i] == "\\newpage ")
        run <- FALSE
      else
        stepfiles <- c(stepfiles, x[i])
      if(i == length(x))
        run <- FALSE
    }
    if(length(id <- grep("\\$f_\\{", stepfiles))) {
      SmoothHyp <- NULL; ok <- FALSE
      for(i in 1:length(id)) {
        term <- paste("sx", collect(stepfiles[id[i]], start = "(", stop = ")"), sep = "")
        nextpart <- strsplit(stepfiles[id[i] + 1], "=")[[1L]]
        if(length(nextpart) == 4L) {
          ok <- TRUE
          lambda <- as.numeric(collect(nextpart[2L], start = " ", stop = " "))
          df <- as.numeric(collect(nextpart[4L], start = " ", stop = " "))
          tmp <- matrix(c(lambda, df), nrow = 1)
          rownames(tmp) <- term
          SmoothHyp <- rbind(SmoothHyp, tmp)
        }
      }
      if(ok)
        colnames(SmoothHyp) <- c("lambda", "df")
      rval$smooth.hyp.step <- SmoothHyp
    }
  }

  return(rval)
}


collect <- function(string, start, stop)
{
  split <- splitme(string)
  take <- NULL; do <- FALSE
  k <- 1
  for(p in split) {
    if(p == start)
      do <- TRUE
    if(do) {
      take <- c(take, p)
    }
    if(k != 1 && p == stop)
      do <- FALSE
    k <- k + 1
  }
  take <- take[take != "$"]

  return(resplit(take))
}


collect2 <- function(string, start, stop)
{
  string <- splitme(string)
  start <- splitme(start)
  stop <- splitme(stop)
  string <- string[string != start & string != stop]

  return(resplit(string))
}

