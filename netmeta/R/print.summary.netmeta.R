print.summary.netmeta <- function(x,
                                  comb.fixed=x$comb.fixed,
                                  comb.random=x$comb.random,
                                  reference.group=x$reference.group,
                                  all.treatments=x$all.treatments,
                                  logscale=FALSE,
                                  header=TRUE,
                                  digits=max(3, .Options$digits - 3),
                                  ...){
  
  
  if (!inherits(x, "summary.netmeta"))
    stop("Argument 'x' must be an object of class \"summary.netmeta\"")
  
  
  k <- x$k
  m <- x$m
  n <- x$n
  sm <- x$sm
  
  sm.lab <- sm
  ##
  if (logscale & meta:::is.relative.effect(sm))
    sm.lab <- paste("log", sm, sep="")

  ci.lab <- paste(round(100*x$level.comb, 1), "%-CI", sep="")

  
  TE.fixed    <- x$fixed$TE
  seTE.fixed  <- x$fixed$seTE
  lowTE.fixed <- x$fixed$lower
  uppTE.fixed <- x$fixed$upper
  ##
  TE.random    <- x$random$TE
  seTE.random  <- x$random$seTE
  lowTE.random <- x$random$lower
  uppTE.random <- x$random$upper
  ##
  if (!is.null(x$seq)){
    TE.fixed <- TE.fixed[x$seq, x$seq]
    seTE.fixed <- seTE.fixed[x$seq, x$seq]
    lowTE.fixed <- lowTE.fixed[x$seq, x$seq]
    uppTE.fixed <- uppTE.fixed[x$seq, x$seq]
    ##
    TE.random <- TE.random[x$seq, x$seq]
    seTE.random <- seTE.random[x$seq, x$seq]
    lowTE.random <- lowTE.random[x$seq, x$seq]
    uppTE.random <- uppTE.random[x$seq, x$seq]
  }
  
  
  noeffect <- 0
  ##
  if (!logscale & meta:::is.relative.effect(sm)){
    noeffect <- 1
    ##
    TE.fixed    <- exp(TE.fixed)
    lowTE.fixed <- exp(lowTE.fixed)
    uppTE.fixed <- exp(uppTE.fixed)
    ##
    TE.random <- exp(TE.random)
    lowTE.random <- exp(lowTE.random)
    uppTE.random <- exp(uppTE.random)
  }
  
  TE.fixed    <- round(TE.fixed, digits)
  lowTE.fixed <- round(lowTE.fixed, digits)
  uppTE.fixed <- round(uppTE.fixed, digits)
  pTE.fixed   <- x$fixed$p
  zTE.fixed   <- round(x$fixed$z, digits)
  ##
  TE.random    <- round(TE.random, digits)
  lowTE.random <- round(lowTE.random, digits)
  uppTE.random <- round(uppTE.random, digits)
  pTE.random   <- x$random$p
  zTE.random   <- round(x$random$z, digits)
  ##
  I2 <- x$I2
  
  
  if (header)
    matitle(x)
  
  if (reference.group!="" & missing(all.treatments))
    all.treatments <- FALSE
  
  if (reference.group !="")
    reference.group <- setref(reference.group, rownames(TE.fixed))
  
  if (comb.fixed|comb.random){
    cat(paste("Number of studies: k=", k, "\n", sep=""))
    cat(paste("Number of treatments: n=", n, "\n", sep=""))
    cat(paste("Number of pairwise comparisons: m=", m, "\n", sep=""))
    
    if (comb.fixed){
      if (all.treatments | reference.group!="")
        cat("\nFixed effect model\n")
      if (all.treatments){
        cat("\nTreatment estimate (sm='", sm.lab, "'):\n", sep="")
        TEf <- format(TE.fixed)
        if (all(diag(TE.fixed)==noeffect))
          diag(TEf) <- "."
        prmatrix(TEf, quote=FALSE, right=TRUE)
        cat("\nLower ", 100*x$fixed$level, "%-confidence limit:\n", sep="")
        lowTEf <- format(lowTE.fixed)
        if (all(diag(lowTE.fixed)==noeffect))
          diag(lowTEf) <- "."
        prmatrix(lowTEf, quote=FALSE, right=TRUE)
        cat("\nUpper ", 100*x$fixed$level, "%-confidence limit:\n", sep="")
        uppTEf <- format(uppTE.fixed)
        if (all(diag(uppTE.fixed)==noeffect))
          diag(uppTEf) <- "."
        prmatrix(uppTEf, quote=FALSE, right=TRUE)
      }
      if (reference.group!=""){
        if (all(colnames(TE.fixed)!=reference.group))
          stop(paste("Argument 'reference.group' must match any of the following values: ",
                     paste(paste("'", colnames(TE.fixed), "'", sep=""),
                           collapse=" - "), sep=""))
        ##
        TE.fixed.b <- TE.fixed[,colnames(TE.fixed)==reference.group]
        lowTE.fixed.b <- lowTE.fixed[,colnames(lowTE.fixed)==reference.group]
        uppTE.fixed.b <- uppTE.fixed[,colnames(uppTE.fixed)==reference.group]
        ##
        res <- cbind(format.TE(TE.fixed.b, na=TRUE),
                     p.ci(format(lowTE.fixed.b), format(uppTE.fixed.b)))
        dimnames(res) <-
          list(colnames(TE.fixed), c(sm.lab, ci.lab))
        ##
        if (TE.fixed.b[rownames(res)==reference.group]==noeffect)
          res[rownames(res)==reference.group,] <- c(".", ".")
        
        cat("\nTreatment estimate (sm='", sm.lab,
            "', reference.group='", reference.group, "'):\n", sep="")
        
        prmatrix(res, quote=FALSE, right=TRUE)
      }
    }
    
    if (comb.random){
      if (all.treatments | reference.group!="")
        cat("\nRandom effects model\n")
      if (all.treatments){
        cat("\nTreatment estimate (sm='", sm.lab, "'):\n", sep="")
        TEr <- format(TE.random)
        if (all(diag(TE.random)==noeffect))
          diag(TEr) <- "."
        prmatrix(TEr, quote=FALSE, right=TRUE)
        cat("\nLower ", 100*x$random$level, "%-confidence limit:\n", sep="")
        lowTEr <- format(lowTE.random)
        if (all(diag(lowTE.random)==noeffect))
          diag(lowTEr) <- "."
        prmatrix(lowTEr, quote=FALSE, right=TRUE)
        cat("\nUpper ", 100*x$random$level, "%-confidence limit:\n", sep="")
        uppTEr <- format(uppTE.random)
        if (all(diag(uppTE.random)==noeffect))
          diag(uppTEr) <- "."
        prmatrix(uppTEr, quote=FALSE, right=TRUE)
      }
      if (reference.group!=""){
        if (all(colnames(TE.random)!=reference.group))
          stop(paste("Argument 'reference.group' must match any of the following values: ",
                     paste(paste("'", colnames(TE.random), "'", sep=""),
                           collapse=" - "), sep=""))
        ##
        TE.random.b <- TE.random[,colnames(TE.random)==reference.group]
        lowTE.random.b <- lowTE.random[,colnames(lowTE.random)==reference.group]
        uppTE.random.b <- uppTE.random[,colnames(uppTE.random)==reference.group]
        ##
        res <- cbind(format.TE(TE.random.b, na=TRUE),
                     p.ci(format(lowTE.random.b), format(uppTE.random.b)))
        dimnames(res) <- list(colnames(TE.fixed), c(sm.lab, ci.lab))
        ##
        if (TE.random.b[rownames(res)==reference.group]==noeffect)
          res[rownames(res)==reference.group,] <- c(".", ".")
        
        
        cat("\nTreatment estimate (sm='", sm.lab,
            "', reference.group='", reference.group, "'):\n", sep="")
        
        prmatrix(res, quote=FALSE, right=TRUE)
      }
    }
    zlab <- "z"
    
    
    if (!is.na(x$tau))
      cat(paste("\nQuantifying heterogeneity/inconsistency:\n",
                if (x$tau^2 > 0 & x$tau^2 < 0.0001)
                paste("tau^2", meta:::format.tau(x$tau^2))
                else
                paste("tau^2 = ",
                      ifelse(x$tau==0,
                             "0",
                             format(round(x$tau^2, 4), 4, nsmall=4, scientific=FALSE)),
                      sep=""),
                paste("; I^2 = ", round(I2, 1), "%",
                      "",
                      ##ifelse(FALSE,
                      ##       p.ci(paste(round(100*lowI2, 1), "%", sep=""),
                      ##            paste(round(100*uppI2, 1), "%", sep="")),
                      ##       ""),
                      sep=""),
                "\n", sep=""))
    

    
    if (m > 1){
      
      Qdata <- cbind(round(x$Q, 2), x$df,
                     ifelse(x$df == 0, "--",
                            meta:::format.p(1-pchisq(x$Q, df=x$df))))
      
      dimnames(Qdata) <- list("", c("Q", "d.f.", "p-value"))
      ##
      cat("\nTest of heterogeneity/inconsistency:\n")
      prmatrix(Qdata, quote=FALSE, right=TRUE, ...)
      ##
      ##
    }
  }
  
  invisible(NULL)
}
