
##NOTE 
# with tsmethod = "stsm", the element "xreg" could be defined in the 
# object "stsm" instead of using argument "xreg" in "tso0". However, 
# it is more convenient to let the function "stsm::stsmFit" handle this element;
# in this way, the same arguments are passed to "stats::arima" and "stsmFit" 
# and the code is simplified here avoiding "if" statements depending on "tsmethod"

tso <- function(y, xreg = NULL, cval = NULL, delta = 0.7, n.start = 50,
  types = c("AO", "LS", "TC"), # c("IO", "AO", "LS", "TC", "SLS")
  maxit = 1, maxit.iloop = 4, cval.reduce = 0.14286, 
  remove.method = c("en-masse", "bottom-up"),
  remove.cval = NULL, 
  tsmethod = c("auto.arima", "arima", "stsm"), 
  args.tsmethod = NULL, args.tsmodel = NULL, logfile = NULL)
{
  tsmethod <- match.arg(tsmethod)
  remove.method <- match.arg(remove.method)
  attr.y <- attributes(y)
  n <- length(y) #attr.y$tsp[2]
  yname <- deparse(substitute(y))
  #stopifnot(is.ts(y))

  if (!is.null(args.tsmethod$xreg))
  {
    if (is.null(xreg))
    {
      # check if external regressors were defined through "args.tsmethod" 
      # instead of using argument "xreg"
      xreg <- args.tsmethod$xreg
      args.tsmethod$xreg <- NULL # this removes element "xreg" from the list
    } else {
      # check if external regressors were defined both 
      # in "xreg" and "args.tsmethod$xreg" but with different values
      if (!identical(xreg, args.tsmethod$xreg))
      {
        warning(paste("non-null", sQuote("args.tsmethod$xreg"), "was ignored; argument ", 
        sQuote("xreg"), "was used instead"))
      } # else # "xreg" was defined twice with the same values (no warning)
      args.tsmethod$xreg <- NULL # this removes element "xreg" from the list
    }
  }

  if (tsmethod == "stsm")
  {
    if (is.null(args.tsmodel$model))
      args.tsmodel$model <- ifelse(frequency(y) == 1, "local-level", "BSM")   

##FIXME these defaults only if stsm.method = "maxlik.fd.scoring"

    if (is.null(args.tsmodel$ssd))
      args.tsmodel$ssd <- TRUE
    if (is.null(args.tsmodel$sgfc))
      args.tsmodel$sgfc <- TRUE
    # let "stsm::stsmFit" handle "xreg", not here
    y <- do.call("stsm.model", args = c(list(y = y), args.tsmodel))
    #ylist <- list(m = m)
  } #else
    #ylist <- list(x = y) # m <- y

  # if "ylist" or "m <- y" were used, then the "if" statement below where "fit" is 
  # created could be avoided using "do.call(tsmethod, args = c(x = m, list())"
  # but this involves storing two identical objects ("y" and "m" or "ylist")

  # default arguments

  if (is.null(args.tsmethod))
  {
    args.tsmethod <- switch(tsmethod,
      "auto.arima" = list(allowdrift = FALSE, ic = "bic"),
      "arima" = list(order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1))),
      "stsm" = list(stsm.method = "maxlik.td.optim", method = "L-BFGS-B",
        KF.version = "KFKSDS", KF.args = list(P0cov = TRUE), gr = "numerical")) #hessian = TRUE
      #list(stsm.method = "maxlik.fd.scoring", step = NULL, information = "expected"))
  }

  # default critical value
  # the same is done in functions "locate.outliers.oloop" and "remove.outliers"
  # "cval" is passed as a non-null value from tso() to those functions
  # but keep there this block so that default value is used when those functions 
  # are called outside tso()

  if (is.null(cval))
  {
    #n <- length(y)
    if (n <= 50) {
      cval <- 3
    } else 
    if (n >= 450) {
      cval <- 4
    } else
      cval <- round(3 + 0.0025 * (n - 50), 2)
  }

  cval0 <- cval
  if (is.null(remove.cval))
    remove.cval <- cval

  # "res0" is used below to generate the output, 
  # "res" is overwritten until no more outliers are found 
  # "res0" is also used if maxit = 1

  res0 <- res <- tso0(x = y, xreg = xreg, cval = cval, 
    delta = delta, n.start = n.start,
    types = types, maxit.iloop = maxit.iloop, 
    remove.method = remove.method, remove.cval = remove.cval,
    tsmethod = tsmethod, args.tsmethod = args.tsmethod, 
    logfile = logfile)

  moall <- res$outliers
  outtimes <- res$times

  iter <- 1
  cval <- round(cval * (1 - cval.reduce), 2)

  while (iter < maxit)
  {
##FIXME see move res0 <- res after if(...) break

if (tsmethod == "stsm")
{
##FIXME TODO create stsm object based on res$yadj as done above
  warning("currently ", sQuote("maxit"), " > 1 is not allowed for ", sQuote("tsmethod=\"stsm\""))
  break
}
    # save "res" to have a copy of the last fitted model, res$fit;
    # if in the current run no outliers are found then 
    # tso0() does not return the fitted model

    res0 <- res

    res <- tso0(x = res$yadj, xreg = xreg, cval = cval, 
      delta = delta, n.start = n.start,
      types = types, maxit.iloop = maxit.iloop, 
      remove.method = remove.method, remove.cval = remove.cval, 
      tsmethod = tsmethod, args.tsmethod = args.tsmethod, 
      logfile = logfile)

    if (nrow(res$outliers) == 0)
      break
    
    moall <- rbind(moall, res$outliers)
    outtimes <- c(outtimes, res$times)
    
    iter <- iter + 1
  }

  # final model given the detected outliers

  if (nrow(moall) > 0)
  {
    pars <- switch(tsmethod, 
      "auto.arima" = , "arima" = coefs2poly(coef(res0$fit), res0$fit$arma, TRUE),
      "stsm" = stsm::char2numeric(res0$fit$model))

    xreg.outl <- outliers.effects(mo = moall, n = n, weights = FALSE, delta = delta, 
      pars = pars, n.start = n.start, freq = frequency(y))
  } else 
    xreg.outl <- NULL

  # all regressors
  # xreg: input regressor variables such as calendar effects (if any)
  # xreg.outl: outliers regressor variables detected above (if any)

  xregall <- cbind(xreg, xreg.outl)
  nms.outl <- colnames(xreg.outl)
  colnames(xregall) <- c(colnames(xreg), nms.outl)

  ##NOTE
  # rerunning "auto.arima" may not be necessary at this point

  if (tsmethod == "stsm") {
    fit <- do.call("stsmFit", args = c(list(x = y, xreg = xregall), args.tsmethod))
  } else {
    fit <- do.call(tsmethod, args = c(list(x = y, xreg = xregall), args.tsmethod))
    # this is for proper printing of results from "auto.arima" and "arima"
    fit$series <- yname
  }

  if (!is.null(xreg.outl))
  {
    id <- colnames(xreg.outl)
    if (tsmethod == "stsm")
    {
##FIXME TODO 
#if xregall!=xreg.outl (i.e. argument xreg is not NULL)
#       xregcoefs <- fit$xreg$coef
#       stde <- fit$xreg$stde
#       if (is.null(stde))
#         stde <- sqrt(diag(vcov(fit, type = "optimHessian")))
      xregcoefs <- fit$pars[id]
      tstats <- xregcoefs / fit$std.errors[id]
    } else { # method "auto.arima", "arima"
      xregcoefs <- coef(fit)[id]
      tstats <- xregcoefs / sqrt(diag(fit$var.coef)[id])
    }

    moall[,"coefhat"] <- xregcoefs
    moall[,"tstat"] <- tstats

    oeff <- xreg.outl %*% cbind(xregcoefs)
    attributes(oeff) <- attr.y #attributes(y)

    yadj <- if(is.ts(y)) y - oeff else y@y - oeff

  } else { # no outliers detected
    oeff <- NULL
    yadj <- if(is.ts(y)) y else y@y
  }

  structure(list(outliers = moall, y = if(is.ts(y)) y else y@y, yadj = yadj, 
    cval = cval0, fit = fit, effects = oeff, times = outtimes), 
    class = "tsoutliers")
}

tso0 <- function(x, xreg = NULL, cval = 3.5, delta = 0.7, n.start = 50,
  types = c("AO", "LS", "TC"), maxit.iloop = 4, 
  remove.method = c("en-masse", "bottom-up"),
  remove.cval = NULL, 
  tsmethod = c("auto.arima", "arima", "stsm"), args.tsmethod = NULL,
  args.tsmodel = NULL, logfile = NULL)
{
  # "x" can be either a "ts" object or a "stsm" object;
  # if !inherits(x, "stsm") then two identical objects are stored ("x" and "y")

  y <- if(is.ts(x)) { x } else x@y

  #remove.method <- match.arg(remove.method)
  #tsmethod <- match.arg(tsmethod)
  #remove.method <- match.arg(remove.method)
  fitmethod <- gsub("stsm", "stsmFit", tsmethod)

  if (is.null(remove.cval))
    remove.cval <- cval

  # fit time series model

  fit <- do.call(fitmethod, args = c(list(x = x, xreg = xreg), args.tsmethod))
  #fit$series <- deparse(substitute(y))

  if (!is.null(logfile))
  {
    cat(paste("model selection:\n"), file = logfile, append = FALSE)
    capture.output(fit, file = logfile, append = TRUE)
  }

  # identify and locate prospective outliers by type
  # given a fitted time series model

  stage1 <- locate.outliers.oloop(y = y, fit = fit, types = types, cval = cval, 
    maxit.iloop = maxit.iloop, delta = delta, n.start = n.start, logfile = logfile)

  # choose and fit the model including the outlier regressors detected so far
  # (the weights of the outliers is fine tuned, to see it 
  # compare 'moall[,"coefhat"]' with 'coef(fit)["oeffi"]') then
  # remove the outliers detected so far if they are not significant in the new model/fit

  if (nrow(stage1$outliers) > 0)
  {
    stage2 <- remove.outliers(x = stage1, y = y, cval = remove.cval, 
      method = remove.method, delta = delta, n.start = n.start, 
      tsmethod.call = fit$call, fdiff = NULL, logfile = logfile)

#moall <- stage2$outliers
stopifnot(ncol(stage2$xreg) == length(stage2$xregcoefs))
  } else 
    stage2 <- list(xreg = NULL, fit = stage1$fit)

  # final outliers and
  # linearized series, original series without outlier effects

  if (!is.null(stage2$xreg))
  {
    # stage2$fit$xreg is not returned by arima()
    moall <- stage2$outliers
    moall[,"coefhat"] <- stage2$xregcoefs
    moall[,"tstat"] <- stage2$xregtstats

    oeff <- stage2$xreg %*% cbind(stage2$xregcoefs)
    attributes(oeff) <- attributes(y)
    yadj <- y - oeff

    moall <- moall[,c("type", "ind", "coefhat", "tstat")]
    outtimes <- time(y)[moall[,"ind"]]
    if (frequency(y) > 1) 
      outseason <- formatC(as.vector(cycle(y)[moall[,"ind"]]), 
        width = 2, flag="0")

    moall <- cbind(moall[,c("type", "ind")], 
      "time" = if (frequency(y) > 1) paste(floor(outtimes), 
        outseason, sep = ":") else outtimes,
      moall[,c("coefhat","tstat")])

    oind <- order(moall[,"ind"])
    moall <- moall[oind,]
    outtimes <- outtimes[oind]
    rownames(moall) <- NULL

  } else { # no outliers detected
    oeff <- NULL
    yadj <- y
    moall <- data.frame(array(dim = c(0, 4)))
    colnames(moall) <- c("type", "ind", "coefhat", "tstat")
    outtimes <- NULL
  }

  if (!is.null(logfile))
  {
    msg <- paste("\nfinal outliers\n")
    cat(msg, file = logfile, append = TRUE)
    capture.output(moall, file = logfile, append = TRUE)
  }

  structure(list(outliers = moall, y = y, yadj = yadj, cval = cval,
    fit = stage2$fit, effects = oeff, times = outtimes), 
    class = "tsoutliers")
}
