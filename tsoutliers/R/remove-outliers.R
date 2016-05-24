
remove.outliers <- function(x, y, cval = NULL, 
  method = c("en-masse", "bottom-up"), 
  delta = 0.7, n.start = 50, tsmethod.call = NULL, 
  fdiff = NULL, logfile = NULL)
{
  if (is.null(tsmethod.call))
    stop(paste(sQuote("tsmethod.call"), "cannot be NULL"))

  # arguments "tsmethod" and "args.tsmethod" (and "args.tsmodel" for "stsm")
  # could be added in order to provide an alternative to argument "tsmethod.call"
  # as way to pass the necessary information;
  # the same that is done in function "tso" could be done here,
  # it would involve dealing again with default values for the fitting method 
  # and for the structural model if "stsm" is selected and the code of this function 
  # would get overly cumbersome; 
  # in addition there is no practical advantage since this function requires 
  # passing in "x" the output of "locate.outliers.oloop" which would have in turn 
  # required calling a function to fit the model and, therefore, the call to that 
  # function should be ready to be passed here

  moall <- x$outliers
  if (is.null(moall) || nrow(moall) == 0)
  {
    #cat("the list of outliers is empty\n")
    return(NULL)
  }

  #tsmethod <- as.character(tsmethod.call[[1]]) 
  #if (tsmethod %in% eval(formals(stsm::stsmFit)$stsm.method))
  #  tsmethod <- "stsm"
  tsmethod <- ifelse(inherits(x$fit$pars, "stsmSS"), 
    "stsm", as.character(tsmethod.call[[1]]))
  method <- match.arg(method)

  # default critical value 
  # (same as in functions "tso" and "locate.outliers.oloop")

  if (is.null(cval))
  {
    n <- length(y)
    if (n <= 50) {
      cval <- 3
    } else 
    if (n >= 450) {
      cval <- 4
    } else
      cval <- round(3 + 0.0025 * (n - 50), 2)
  }

  # regressor variables (other than outliers)
  # "xregfixed" will be either NULL or a matrix with column names

  xregfixed <- tsmethod.call$xreg
  xregfixed.nms <- colnames(xregfixed)
  tsmethod.call$xreg <- NULL

  # outlier regressor variables

  ##NOTE
  #"pars" is used only if IO is considered

  xreg <- outliers.effects(mo = moall, n = x$fit$n, weights = FALSE,
    delta = delta, pars = x$fit$pars, n.start = n.start, freq = frequency(y))

  iter <- 0

  if (method == "en-masse")
  {
    while (TRUE)
    {
      #if (!is.null(tsmethod.call))
      xregall <- cbind(xregfixed, xreg)
      colnames(xregall) <- c(xregfixed.nms, colnames(xreg))

      ##NOTE
      #"eval(as.call(c(as.list(tsmethod.call), list(xreg = xregall))))" 
      # could be used for "auto.arima" as well but the "if" statement below 
      # is required because of the following behaviour
      #
      # a <- arima(x = y, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1)), 
      #   include.mean = FALSE)
      # a$call
      # arima(x = y, order = c(0, 1, 1), seasonal = list(order = c(0, 
      #     1, 1)), include.mean = FALSE)
      # a <- do.call("arima", args = c(list(x = y, order = c(0, 1, 1), 
      #   seasonal = list(order = c(0, 1, 1)), include.mean = FALSE)))
      # class(a$call$x)
      # [1] "ts"
      # identical(a, eval(a$call))
      # [1] TRUE
      #
      # b <- auto.arima(x = y, allowdrift = FALSE, ic = "bic")
      # class(b$call$x)
      # [1] "data.frame"
      # identical(b, eval(b$call))
      # [1] FALSE

      if (tsmethod.call[[1]] == "auto.arima") {
        tsmethod.call$x <- NULL # this could be done outside this loop
        fit <- do.call("auto.arima", args = c(list(x = y), 
          as.list(tsmethod.call[-1]), list(xreg = xregall)))
      } else {
        fit <- eval(as.call(c(as.list(tsmethod.call), list(xreg = xregall))))
      }

      if (!is.null(logfile))
      {
        msg <- paste("\nchoose model and remove outliers, iter:", iter, "\n")
        cat(msg, file = logfile, append = TRUE)
        capture.output(fit, file = logfile, append = TRUE)
      }

      if (tsmethod == "stsm")
      {
        #print(fit, vcov.type = "optimHessian")
        #xregcoefs <- fit$xreg$coef
        #stde <- fit$xreg$stde
        #if (is.null(stde))
        #  stde <- sqrt(diag(vcov(fit, type = "optimHessian")))
        #tstats <- xregcoefs / stde[names(xregcoefs)]
        id <- colnames(xreg)
        xregcoefs <- fit$pars[id]
        tstats <- xregcoefs / fit$std.errors[id]
      } else { # method "auto.arima", "arima"        
        # "xreg" is not returned by arima (fit$call$xreg could be used)
        #id <- colnames(fit$xreg)
        id <- colnames(xreg)
        xregcoefs <- coef(fit)[id]
        tstats <- xregcoefs / sqrt(diag(fit$var.coef)[id])
      }

      # remove outliers if they are not significant

      ref <- which(abs(tstats) < cval)

      if (length(ref) > 0)
      {
        moall <- moall[-ref,]
        # data.matrix() does not keep the column name when ncol=1
        # xreg <- data.matrix(xreg[,-ref])
        xreg <- xreg[,-ref]
        if (is.null(dim(xreg))) {
          xreg <- matrix(xreg, ncol = 1)
          colnames(xreg) <- paste(as.character(moall[,"type"]), moall[,"ind"], sep = "")
        }

        if (!is.null(logfile))
        {
          cat(paste("outliers, iter:", iter, "\n"), file = logfile, append = TRUE)
          capture.output(moall, file = logfile, append = TRUE)
        }
      } else 
        break

      if (nrow(moall) == 0)
        break

      iter  <- iter + 1
    }
  } else
  if (method == "bottom-up")
  {
    ref <- order(abs(moall[,"tstat"]), decreasing = TRUE)
    id <- idrm <- NULL
    xregaux <- matrix(nrow = nrow(xreg), ncol = 0)

    for (i in ref)
    {
      id <- c(id, i)
      if (length(id) == 1) {
        xregaux <- matrix(xreg[,i], ncol = 1, dimnames = list(NULL, colnames(xreg)[i]))
      } else {
stopifnot(length(id) > 1)
        xregaux <- xreg[,id]
      }

      if (!is.null(xregfixed))
      {
        xregall <- cbind(xregfixed, xregaux)
        colnames(xregall) <- c(xregfixed.nms, colnames(xregaux))
      } else xregall <- xregaux

      if (tsmethod.call[[1]] == "auto.arima") {
        tsmethod.call$x <- NULL # this could be done before the current loop
        fit <- do.call("auto.arima", args = c(list(x = y), 
          as.list(tsmethod.call[-1]), list(xreg = xregall)))  #xregaux
      } else {
        fit <- try(eval(as.call(c(as.list(tsmethod.call), list(xreg = xregall)))), #xregaux
          silent = TRUE)
      }

      if (!is.null(logfile))
      {
        msg <- paste("\nchoose model and remove outliers, iter:", iter, "\n")
        cat(msg, file = logfile, append = TRUE)
        capture.output(fit, file = logfile, append = TRUE)
      }

      if (tsmethod == "stsm")
      {
        #print(fit, vcov.type = "optimHessian")
        #xregcoefs <- fit$xreg$coef
        #stde <- fit$xreg$stde
        #if (is.null(stde))
        #  stde <- sqrt(diag(vcov(fit, type = "optimHessian")))
        #tstats <- xregcoefs / stde[names(xregcoefs)]
        nms <- colnames(xregaux)
        xregcoefs <- fit$pars[nms]
        tstats <- xregcoefs / fit$std.errors[nms]
      } else { # method "auto.arima", "arima"
      if (!inherits(fit, "try-error"))
        {
          nms <- colnames(xregaux)
          xregcoefs <- coef(fit)[nms]
          tstats <- xregcoefs / sqrt(diag(fit$var.coef)[nms])
        } else {
          xregcoefs <- NULL
          tstats <- 0
        }
      }

      # if adding a regressor, say "xregi", makes any of the variables in 
      # "xregaux" not significant then "xregi" is removed

      if (any(abs(tstats) < cval))
      {
        id <- id[-length(id)]
        idrm <- c(idrm, i)

        if (i == tail(ref, 1)) {
          xregcoefs <- xregcoefs[-length(xregcoefs)]
          tstats <- tstats[-length(tstats)]
          recover.fit <- TRUE
          fit0 <- fit
        } else
          recover.fit <- FALSE
      } else {
        fit0 <- fit
        recover.fit <- FALSE
      }
    }

    # recover the last fit where the regressors where significant
    # so that the output that is returned in "fit" is in agreement with "moall"

    if (recover.fit)
      fit <- fit0

    if (length(idrm) > 0) {
      #rm(xregaux) #"xregaux" will not be in general a big object
      xreg <- xreg[,-idrm]
      moall <- moall[-idrm,]
    }
  } #else # not necessary, it would be caught by match.arg(method) above
    #stop("unkown method")

  if (nrow(moall) == 0)
  {
stopifnot(ncol(xreg) == 0)
    fit <- xreg <- xregcoefs <- tstats <- NULL
    #xreg <- NULL
  }

  list(xreg = xreg, xregcoefs = xregcoefs, xregtstats = tstats, 
    iter = iter, fit = fit, outliers = moall)
}
