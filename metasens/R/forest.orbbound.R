forest.orbbound <- function(x,
                            comb.fixed=x$x$comb.fixed,
                            comb.random=x$x$comb.random,
                            text.fixed="FE model",
                            text.random="RE model",
                            smlab=NULL,
                            leftcols=c("studlab", "maxbias"),
                            leftlabs=c("Missing\nstudies", "Maximum\nbias"),
                            backtransf=x$backtransf,
                            digits=max(3, .Options$digits - 3),
                            ...){
  
  meta:::chkclass(x, "orbbound")
  
  k <- x$x$k
  sm <- x$x$sm
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.summary.limitmeta"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg="backtransf")
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  
  
  if (length(comb.fixed)==0)
    comb.fixed <- TRUE
  ##
  if (length(comb.random)==0)
    comb.random <- TRUE
  
  
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
  
  
  ci.lab <- paste(round(100*x$x$level.comb, 1), "%-CI", sep="")
  
  
  TE.fixed   <- x$fixed$TE
  seTE.fixed <- rep(x$fixed$seTE, length(TE.fixed))
  ##
  TE.random   <- x$random$TE
  seTE.random <- rep(x$random$seTE, length(TE.random))
  
  
  if (comb.fixed & comb.random){
    TE <- c(TE.fixed, TE.random)
    seTE <- c(seTE.fixed, seTE.random)
    FEvsRE <- c(rep(text.fixed, length(TE.fixed)),
                rep(text.random, length(TE.random)))
  }
  if (comb.fixed & !comb.random){
    TE <- TE.fixed
    seTE <- seTE.fixed
    if (is.null(smlab))
      smlab <- "Fixed effect model"
  }
  if (!comb.fixed & comb.random){
    TE <- TE.random
    seTE <- seTE.random
    if (is.null(smlab))
      smlab <- "Random effects model"
  }
  if (!comb.fixed & !comb.random){
    warning("No forest plot generated as both arguments 'comb.fixed' and 'comb.random' are FALSE")
    return(invisible(NULL))
  }
  
  
  if (comb.fixed & comb.random)
    m1 <- metagen(TE, seTE, sm=sm.lab,
                  byvar=FEvsRE, print.byvar=FALSE,
                  warn=FALSE)
  else
    m1 <- metagen(TE, seTE, sm=sm.lab, warn=FALSE)
  ##
  if (comb.fixed & comb.random){
    m1$studlab <- c(x$k.suspect, x$k.suspect)
    m1$maxbias <- c(x$maxbias, x$maxbias)
    m1$npft.ma <- c(1/mean(1/x$x$n), 1/mean(1/x$x$n))
  }
  else{
    m1$studlab <- x$k.suspect
    m1$maxbias <- x$maxbias
    m1$npft.ma <- 1/mean(1/x$x$n)
  }
  
  
  if (backtransf)
    m1$maxbias <- meta:::backtransf(m1$maxbias, sm, "mean",
                                    m1$npft.ma, warn=FALSE)
  ##
  m1$maxbias <- format(round(m1$maxbias, digits))
  
  forest(m1,
         comb.fixed=FALSE, comb.random=FALSE,
         hetstat=FALSE,
         leftcols=leftcols,
         leftlabs=leftlabs,
         smlab=smlab,
         just.studlab="center",
         weight="same",
         ...)
  
  
  invisible(NULL)
}
