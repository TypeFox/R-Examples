print.limitmeta <- function(x,
                            sortvar,
                            backtransf=x$backtransf,
                            digits=max(3, .Options$digits - 3),
                            header=TRUE, ...){
  
  meta:::chkclass(x, "limitmeta")
  
  
  format.TE <- function(TE, na=FALSE){
    TE <- meta:::rmSpace(TE)
    if (na) res <- format(TE)
    else res <- ifelse(is.na(TE), "", format(TE))
    res
  }
  
  
  sm <- x$sm
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.limitmeta"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg="backtransf")
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  
  
  sm.lab <- sm
  ##
  if (backtransf){
    if (sm=="ZCOR")
      sm.lab <- "COR"
    if (sm %in% c("PFT", "PAS", "PRAW", "PLOGIT", "PLN"))
      sm.lab <- "proportion"
  }
  else 
    if (meta:::is.relative.effect(sm))
      sm.lab <- paste("log", sm, sep="")
  
  
  ci.lab <- paste(round(100*x$level, 1), "%-CI", sep="")
  
  
  TE <- x$TE
  seTE <- x$seTE
  ##
  TE.limit <- x$TE.limit
  seTE.limit <- x$seTE.limit
  
  
  k.all <- length(TE)
  ##
  if (missing(sortvar)) sortvar <- 1:k.all
  ##
  if (length(sortvar) != k.all)
    stop("Arguments 'x' and 'sortvar' have different length")
  
  
  ci.TE <- ci(TE, seTE, level=x$level)
  lowTE <- ci.TE$lower
  uppTE <- ci.TE$upper
  ##
  ci.limit <- ci(TE.limit, seTE.limit, level=x$level)
  lowTE.limit <- ci.limit$lower
  uppTE.limit <- ci.limit$upper
  
  
  if (backtransf){
    ##
    npft.ma <- 1/mean(1/x$x$n)
    ##
    TE    <- meta:::backtransf(TE, sm, "mean",
                               npft.ma, warn=TRUE)
    lowTE <- meta:::backtransf(lowTE, sm, "lower",
                               npft.ma, warn=TRUE)
    uppTE <- meta:::backtransf(uppTE, sm, "upper",
                               npft.ma, warn=TRUE)
    ##
    TE.limit <- meta:::backtransf(TE.limit, sm, "mean",
                                  npft.ma, warn=TRUE)
    lowTE.limit <- meta:::backtransf(lowTE.limit, sm, "lower",
                                     npft.ma, warn=TRUE)
    uppTE.limit <- meta:::backtransf(uppTE.limit, sm, "upper",
                                     npft.ma, warn=TRUE)
  }
  ##
  TE    <- round(TE, digits)
  lowTE <- round(lowTE, digits)
  uppTE <- round(uppTE, digits)
  ##
  TE.limit    <- round(TE.limit, digits)
  lowTE.limit <- round(lowTE.limit, digits)
  uppTE.limit <- round(uppTE.limit, digits)
  ##
  TE    <- format.TE(TE, na=TRUE)
  lowTE <- format(lowTE, trim=TRUE)
  uppTE <- format(uppTE, trim=TRUE)
  ##
  TE.limit    <- format.TE(TE.limit, na=TRUE)
  lowTE.limit <- format(lowTE.limit, trim=TRUE)
  uppTE.limit <- format(uppTE.limit, trim=TRUE)
  
  
  res <- cbind(" ",
               TE,
               meta:::p.ci(lowTE, uppTE),
               "  ",
               TE.limit,
               meta:::p.ci(lowTE.limit, uppTE.limit))
  ##
  dimnames(res) <-
    list(x$studlab, c("", sm.lab, ci.lab, "", sm.lab, ci.lab))
  ##
  cat("Results for individual studies (left: original data; right: shrunken estimates)\n\n")
  prmatrix(res[order(sortvar),], quote=FALSE, right=TRUE)
  
  
  cat("\n")
  
  
  print(summary(x),
        digits=digits, header=FALSE, backtransf=backtransf)
  
  
  invisible(NULL)
}
