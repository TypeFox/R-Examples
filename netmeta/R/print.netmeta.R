print.netmeta <- function(x,
                          sortvar,
                          level=x$level, level.comb=x$level.comb,
                          comb.fixed=x$comb.fixed, comb.random=x$comb.random,
                          reference.group=x$reference.group,
                          all.treatments=x$all.treatments,
                          details=TRUE, ma=TRUE, logscale=FALSE,
                          digits=max(4, .Options$digits - 3),
                          ...
                          ){
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")
  
  
  k.all <- length(x$TE)
  ##
  if (missing(sortvar)) sortvar <- 1:k.all
  ##
  if (length(sortvar) != k.all)
    stop("'x' and 'sortvar' have different length")
  ##
  ci.lab <- paste(round(100*level, 1), "%-CI", sep="")

  sm <- x$sm
  
  sm.lab <- sm
  ##
  if (logscale & meta:::is.relative.effect(sm))
    sm.lab <- paste("log", sm, sep="")

  
  matitle(x)
  
  
  if (details){

    cat(paste("Original data",
              ifelse(any(x$narms>2),
                     " (with adjusted standard errors for multi-arm studies)",
                     ""),
              ":\n\n", sep=""))
    
    res <- data.frame(treat1=x$treat1,
                      treat2=x$treat2,
                      TE=format.TE(round(x$TE, digits), na=TRUE),
                      seTE=format(round(x$seTE, digits)))
    ##
    if (any(x$narms>2))
      res$seTE.adj <- format(round(x$seTE.adj, digits))
    ##
    res$studlab <- x$studlab
    ##
    if (any(x$narms>2)){
      tdata1 <- data.frame(studlab=as.character(x$studies),
                           narms=x$narms)
      res$OrDeR <- 1:dim(res)[[1]]
      res <- merge(res, tdata1,
                   by="studlab", all.x=TRUE, all.y=FALSE,
                   sort=FALSE)
      res <- res[order(res$OrDeR),]
      res$multiarm <- ifelse(res$narms>2, "*", "")
      res$OrDeR <- NULL
      res$studlab <- NULL
      res <- as.matrix(res)
    }
    else
      res$studlab <- NULL
    ##
    dimnames(res)[[1]] <- x$studlab
    
    prmatrix(res[order(sortvar),], quote=FALSE, right=TRUE)
    cat("\n")
    
    cat("Number of treatment arms (by study):\n")
    prmatrix(data.frame(narms=x$narms, row.names=x$studies),
             quote=FALSE, right=TRUE)
    cat("\n")
  }
  
  
  tsum <- summary(x, level=level, level.comb=level.comb, warn=FALSE)
  ##
  TE.f    <- tsum$comparison.nma.fixed$TE
  lowTE.f <- tsum$comparison.nma.fixed$lower
  uppTE.f <- tsum$comparison.nma.fixed$upper
  ##
  if (!logscale & meta:::is.relative.effect(sm)){
    TE.f    <- exp(TE.f)
    lowTE.f <- exp(lowTE.f)
    uppTE.f <- exp(uppTE.f)
  }
  ##
  TE.f <- round(TE.f, digits)
  lowTE.f <- round(lowTE.f, digits)
  uppTE.f <- round(uppTE.f, digits)
  ##
  ##
  TE.r    <- tsum$comparison.nma.random$TE
  lowTE.r <- tsum$comparison.nma.random$lower
  uppTE.r <- tsum$comparison.nma.random$upper
  ##
  if (!logscale & meta:::is.relative.effect(sm)){
    TE.r    <- exp(TE.r)
    lowTE.r <- exp(lowTE.r)
    uppTE.r <- exp(uppTE.r)
  }
  ##
  TE.r <- round(TE.r, digits)
  lowTE.r <- round(lowTE.r, digits)
  uppTE.r <- round(uppTE.r, digits)
  
  
  res.f <- cbind(x$treat1, x$treat2,
                 format.TE(TE.f, na=TRUE),
                 p.ci(format(lowTE.f), format(uppTE.f)),
                 if (comb.fixed) format(round(x$Q.fixed, 2)),
                 if (comb.fixed) format(round(x$leverage.fixed, 2)))
  dimnames(res.f) <-
    list(x$studlab, c("treat1", "treat2",
                      sm.lab, ci.lab,
                      if (comb.fixed) "Q",
                      if (comb.fixed) "leverage"))
  
  
  res.r <- cbind(x$treat1, x$treat2,
                 format.TE(TE.r, na=TRUE),
                 p.ci(format(lowTE.r), format(uppTE.r)))
  dimnames(res.r) <-
    list(x$studlab, c("treat1", "treat2", sm.lab, ci.lab))
  
  
  if (comb.fixed){
    cat("Results (fixed effect model):\n\n")
    
    prmatrix(res.f[order(sortvar),], quote=FALSE, right=TRUE)
    
    cat("\n")
  }
  
  if (comb.random){
    cat("Results (random effects model):\n\n")
    
    prmatrix(res.r[order(sortvar),], quote=FALSE, right=TRUE)
    
    cat("\n")
  }
  
  if (reference.group!="" & missing(all.treatments))
    all.treatments <- FALSE
  
  
  if (reference.group !="")
    reference.group <- setref(reference.group, rownames(x$A.matrix))
  
  
  if (ma)
    print(tsum, digits=digits,
          comb.fixed=comb.fixed, comb.random=comb.random,
          logscale=logscale,
          all.treatments=all.treatments,
          reference.group=reference.group,
          header=FALSE)
  
  invisible(NULL)
}
