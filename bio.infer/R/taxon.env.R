"taxon.env" <-
  function(form, bcnt, envdata,  bcnt.siteid, bcnt.abndid,
                      env.siteid, tlevs = "all", dumpdata = FALSE) {

  if (!(bcnt.siteid %in% names(bcnt))) {
    stop(bcnt.siteid, " not found in bcnt")
  }
  if (!(bcnt.abndid %in% names(bcnt))) {
    stop(bcnt.abndid, " not found in bcnt")
  }
  if (!(env.siteid %in% names(envdata))) {
    stop(env.siteid, " not found in envdata")
  }
  
  # Reformat form as a regression formula
  formstr <- as.character(form)[2]
  split0 <- strsplit(formstr, "\\+")[[1]]
  for (i in 1:length(split0)) {
    split0[i] <- sub(" +$", "", split0[i])
    split0[i] <- sub("^ +", "", split0[i])
  }

  w <- regexpr("\\^", split0)
  w2 <- regexpr("\\*", split0)
  xvar <- character(0)
  for (i in 1:length(split0)) {
    if (w[i] != -1 ) {
      xvar <- c(xvar, substring(split0[i], 1, w[i]-1))
      split0[i] <- paste("I(", split0[i],")", sep = "")
    }
    else {
      if (w2[i] != -1) {
        xvar <- c(xvar, substring(split0[i], 1, w2[i]-1))
        xvar <- c(xvar, substring(split0[i], w2[i]+1, nchar(split0[i])))
      }
      else {
        xvar <- c(xvar, split0[i])
      }
    }
  }

  for (i in 1:length(xvar)) {
    xvar[i] <- sub(" +$", "", xvar[i])
    xvar[i] <- sub("^ +", "", xvar[i])
  }

  xvar <- unique(xvar)

  for (i in 1:length(xvar)) {
    if (!(xvar[i] %in% names(envdata))) {
      stop("Variable ",xvar[i], " not found in environmental data.")
    }
  }
  dfenv <- na.omit(envdata[, c(env.siteid, xvar)])

  # A partial check of the match between the variable list and the formula...
  varinc <- rep(FALSE, times = length(split0))
  for (i in 1:length(xvar)) {
    w1 <- regexpr(xvar[i], split0)
    varinc <- varinc | (w1 != -1)
  }
  if (sum(varinc) != length(split0)) {
    stop("Formula is not consistent with variable list.")
  }
  
  form0 <- paste(split0, collapse = "+")
  form0 <- paste("resp ~", form0, sep = "")

  cat("Model formula: ", form0, "\n")

  # set cutoff number of observations as 10 times the number of df
  cutoff <- (length(split0)+1)*10

  cat("Minimum number of occurrences: ", cutoff, "\n")

  if (tlevs[1] == "all") {
    tlevs.sel <- names(bcnt)[4:length(names(bcnt))]
  }
  else {
    for (i in 1:length(tlevs)) {
      if (! (tlevs[i] %in% names(bcnt))) {
        stop(tlevs[i], "not found in bcnt")
      }
    }
    tlevs.sel <- tlevs
  }

  bcnt0 <- merge(bcnt, dfenv, by.x = bcnt.siteid, by.y = env.siteid)

  xlims <- as.list(rep(NA, times = length(xvar)))
  for (i in 1:length(xvar)) {
    r0 <- range(bcnt0[, xvar[i]])
    xlims[[i]] <- r0
    dfenv[,xvar[i]] <- (dfenv[,xvar[i]] - r0[1])/diff(r0)
  }

  if (is.factor(bcnt0[,bcnt.siteid])) {
    sitenames.u <- sort(unique(levels(bcnt0[,bcnt.siteid])[bcnt0[,bcnt.siteid]]))
  }
  else{
    sitenames.u <- sort(unique(bcnt0[,bcnt.siteid]))
  }

  numid <- is.numeric(sitenames.u)
    
  nsamp <- length(sitenames.u)
  tnames.sav <- character(0)
  roc <- numeric(0)
  luniq <- function(x) length(unique(x))

  tcount <- 0
  nummod <- rep(NA, length(tlevs.sel))

  raw.data <- as.list(rep(NA, times = length(tlevs.sel)))

  for (i in 1:length(tlevs.sel)) {
    numocc <- tapply(bcnt0[, bcnt.siteid], bcnt0[,tlevs.sel[i]], luniq)
    numocc2 <- nsamp - numocc
    numocc3 <- pmin(numocc, numocc2)
    incvec <- numocc3 >= cutoff

    tnames.loc <- names(numocc3)[incvec]
    nummod[i] <- length(tnames.loc)

    if (length(tnames.loc) > 0) {
      ntaxa <- length(tnames.loc)
    # Build site species matrix for only taxa with enough observations
      ntaxa <- length(tnames.loc)
      abund.all <- bcnt0[, bcnt.abndid]
      ss1 <- matrix(0, ncol = ntaxa, nrow = nsamp)
      for (j in 1:ntaxa) {
        incvec <- bcnt0[, tlevs.sel[i]] == tnames.loc[j]
        incvec[is.na(incvec)] <- FALSE
        sitenames.g <- bcnt0[incvec,bcnt.siteid]
        abund.g <- bcnt0[incvec,bcnt.abndid]
        abund.s <- tapply(abund.g, sitenames.g, sum)
        sitenames.s <- names(abund.s)
        if (numid) {
          sitenames.s <- as.numeric(sitenames.s)
        }
        for (k in 1:length(sitenames.s)) {
          ss1[match(sitenames.s[k], sitenames.u), j] <- abund.s[k]
        }
      }
      ss <- data.frame(sitenames.u, ss1)
      names(ss) <- c("SITEID", tnames.loc)
      df1 <- merge(ss, dfenv, by.x = "SITEID", by.y = env.siteid)

      raw.data[[i]] <- df1

      for (j in 1:ntaxa) {
        resp <- df1[, tnames.loc[j]] > 0
        resp[is.na(resp)] <- FALSE
        options(warning.expression = expression())
        withCallingHandlers(mod1 <- glm(as.formula(form0),
                                        data = df1, family = "binomial"),
                            warning = function(w) cat("Warning:", conditionMessage(w), "\n\tfor", tnames.loc[j], "\n")
                            )
        options(warning.expression = NULL)
        tcount <- tcount + 1

        predout <- predict(mod1, type = "response")

        x <- predout[resp]
        y <- predout[! resp]

        rocmat <- matrix(NA, nrow = length(x), ncol = length(y))
        for (j in 1:length(x)) {
          rocmat[j,] <- as.numeric(x[j] > y)
        }

        roc <- c(roc, round(sum(rocmat)/(length(x)*length(y)), digits=3))
        
        if (tcount == 1) {
          csave <- as.vector(coef(mod1))
        }
        else{
          csave <- rbind(csave, as.vector(coef(mod1)))
        }
      }
      tnames.sav <- c(tnames.sav, tnames.loc)
    }
  }

  dfreport <- data.frame(tlevs.sel, nummod)
  names(dfreport) <- c("TAXON.LEVEL", "NUM.MODS")
  cat("Number of taxa modeled:\n")
  print(dfreport)
  

  if (dumpdata) {
    coef0 <- list(tnames = tnames.sav, csave = csave, xvar= xvar,
                  xlims = xlims, form = form, roc = roc,
                  raw.data = raw.data)
  }
  else { 
    coef0 <- list(tnames = tnames.sav, csave = csave, xvar= xvar,
                  xlims = xlims,form = form, roc = roc)
  }

  return(coef0)
  
}


