
ch.test <- function(x, type = c("dummy", "trigonometric"), 
  lag1 = FALSE, NW.order = NULL, sid = NULL, xreg = NULL, 
  pvalue = c("RS", "raw"), rs.nobsreg = 13)
{
  ch.test0 <- function(id)
  {
    # robust covariance matrix estimator 
    # scaled by the seasonal components
    # based on function "sandwich::NeweyWest"
    # (other options used in package "sandwich" could be considered here)
    # OmegafHat <- crossprod(SD * ehat) / n
    SDe <- SD[,id,drop=FALSE] * ehat
    if (NW.order > 0)
    {
      weights <- 1 - seq_len(NW.order) / (NW.order + 1)
      OmegafHat <- 0.5 * crossprod(SDe)
      for (i in seq_along(weights))
        OmegafHat <- OmegafHat + 
        weights[i] * crossprod(SDe[seq_len(n-i),,drop=FALSE], SDe[-seq_len(i),,drop=FALSE])
      OmegafHat <- OmegafHat + t(OmegafHat)
    } else
      OmegafHat <- crossprod(SDe)
    OmegafHat <- OmegafHat / n
    # auxiliary elements
    # cumulative sum by columns in "SDe" ("SD" multiplied by "ehat" by columns)
    Fhat <- apply(SDe, MARGIN = 2, FUN = cumsum)
    # summation of crossproducts of "Fhat"
    Fhat.cp <- matrix(0, nrow = nrow(OmegafHat), ncol = ncol(OmegafHat))
    for (i in seq_len(n))
      Fhat.cp <- Fhat.cp + tcrossprod(Fhat[i,])
    # test statistic
    sum(diag(chol2inv(chol(OmegafHat)) %*% Fhat.cp)) / n^2
  }

  data.name <- deparse(substitute(x))
  type <- match.arg(type)
  pvalue <- match.arg(pvalue)
  isNullxreg <- is.null(xreg)

  n <- length(x)
  if (!isNullxreg && NROW(xreg) != n)
    stop("wrong dimension of argument ", sQuote("xreg"))
  S <- frequency(x)
  if (S < 2)
    stop("the time series is not seasonal")

  if (is.null(NW.order))
    NW.order <- round(S * (n/100)^0.25)

  # indicator variable for the target seasonal dummies or cycles

  if (is.null(sid)) {
    sid <- "all"
  } else {
    #if (!identical(sid, "all"))
    if (is.numeric(sid))
    {
      if (type == "trigonometric")
      {
        if (length(sid) != floor(S/2))
          stop("wrong length of argument ", sQuote("sid"))
        tmp <- head(sid, -1)
        id <- which(c(rbind(tmp, tmp)) == 1)
        if (tail(sid, 1) == 1)
          id <- c(id, S-1)
      } else {
        if (any(!(sid %in% seq_len(S)))) # assumed that no duplicates are defined
          stop("wrong definition of argument ", sQuote("sid"))
        id <- sid
      }
    } else
      if (!(sid %in% c("all", "joint")))
        stop("wrong definition of argument ", sQuote("sid"))
  }

  # create target regressor variables (seasonal dummies or seasonal cycles)

  switch(type,
    "dummy" = { # seasonal dummies
      #SD <- seasonal.dummies(x)
      SD <- do.call("rbind", replicate(ceiling(n/S), diag(S), simplify = FALSE))
      SD <- ts(SD, frequency = S, start = c(start(x)[1], 1))
      # ignore warning "'end' value not changed"
      SD <- suppressWarnings(window(SD, start = start(x), end = end(x)))  
    },
    "trigonometric" = { # seasonal cycles
      Sh <- floor(S/2)
      isSeven <- as.numeric(S %% 2 == 0)
      if (S %in% c(4, 12) && n/S <= 50) {
        SD <- .SDtrig[[as.character(S)]][seq_len(n),]
      } else {
        #SD <- seasonal.cycles(x)
        tmp <- matrix(seq_len(n), nrow = Sh-isSeven, ncol = n, byrow = TRUE)
        seqsm1 <- seq_len(nrow(tmp))
        tmp <- (2 * seqsm1 * pi / S) * tmp
        SD <- rbind(cos(tmp), sin(tmp))
        SD <- t(SD[c(rbind(seqsm1, seqsm1 + nrow(tmp))),])
        if ((S %% 2) == 0)
          SD <- cbind(SD, rep(c(-1, 1), len = n))
        #SD <- ts(SD, frequency = S, start = c(start(x)[1], 1))
      }
    }
  ) # switch
  colnames(SD) <- paste0("SD", seq_len(ncol(SD)))

  # arrange other possible exogenous variables

  if (lag1) {
    SD <- SD[-1,] # updating SD is required as well as it will be used by "ch.test0"
    xreg <- cbind(lag1 = x[-n], SD, xreg[-1,])
    x <- x[-1]
    n <- n - 1
  } else
    xreg <- cbind(SD, xreg)

  ##NOTE
  # in principle it is not a good idea to define the intercept in "xreg" and 
  # use lm(x ~ 0 + xreg) because stats::summary.lm uses attr(z$terms, "intercept")
  # to compute the R-squared, but here the R-squared is not used

  if (type == "trigonometric")
    xreg <- cbind(c = 1, xreg)

  # fit regression model and get residuals

  fit <- lm(x ~ 0 + xreg)
  ehat <- residuals(fit)

  # used with pvalue = "RS"
  Nc <- n - ncol(xreg)

  # test statistics

  if (identical(sid, "all"))
  {
    if (type == "dummy")
    {
      stat <- matrix(nrow = S+1, ncol = 2)
      id <- 0
      for (i in seq_len(S))
      {
        id <- id + 1
        stat[i,1] <- ch.test0(id)
        stat[i,2] <- switch(pvalue,
          #"raw" = 1 - .CH.cvals[[1]](stat[i,1]),
          "raw" = uroot.raw.pvalue(stat[i,1], "CH", 1),
          "RS" = ch.rs.pvalue(stat[i,1], "dummy", lag1, S, Nc, rs.nobsreg, 1))
      }
      stat[S+1,1] <- ch.test0(seq_len(S))
      ##NOTE to mention in documentation
      #"RS" p-value is not available for the joint test with dummies
      stat[S+1,2] <- uroot.raw.pvalue(stat[S+1,1], "CH", S)
      colnames(stat) <- c("statistic", "p-value")
      rownames(stat) <- c(switch(as.character(S),
        "4" = paste0("Quarter", seq_len(S)), "12" = month.abb, 
        paste0("Season", seq_len(S))), "joint")

    } else { #type == "trigonometric"
      Shp1 <- Sh + 1
      stat <- matrix(nrow = Shp1, ncol = 2)
      id <- c(-1, 0)
      for (i in seq_len(Sh-isSeven))
      {
        id <- id + 2
        stat[i,1] <- ch.test0(id)
        stat[i,2] <- switch(pvalue,
          "raw" = uroot.raw.pvalue(stat[i,1], "CH", 2),
          "RS" = ch.rs.pvalue(stat[i,1], "trigonometric", lag1, S, Nc, rs.nobsreg, 2))
      }
      if (isSeven)
      {
        stat[Sh,1] <- ch.test0(S-1)
        stat[Sh,2] <- switch(pvalue,
          "raw" = uroot.raw.pvalue(stat[Sh,1], "CH", 1),
          "RS" = ch.rs.pvalue(stat[Sh,1], "trigonometric", lag1, S, Nc, rs.nobsreg, 1))
      }
      stat[Shp1,1] <- ch.test0(seq_len(S-1))
      stat[Shp1,2] <- switch(pvalue,
        "raw" = uroot.raw.pvalue(stat[Shp1,1], "CH", S-1),
        "RS" = ch.rs.pvalue(stat[Shp1,1], "trigonometric", lag1, S, Nc, rs.nobsreg, S-1))
      colnames(stat) <- c("statistic", "p-value")
      if (isSeven) {
        rownames(stat) <- c(paste0("pi/", Sh), 
          if(Sh > 2) paste0(seq.int(2, Sh-1), "pi/", Sh), "pi", "joint")
      } else 
        rownames(stat) <- c(paste0(seq.int(2, S, 2), "pi/", S), "joint")
    }
  } else # sid != "all"
  if (identical(sid, "joint")) {
    if (type == "dummy")
    {
      stat <- ch.test0(seq_len(S))
      stat <- cbind("joint" = stat, "p-value" = uroot.raw.pvalue(stat, "CH", S))
    } else { # type trigonometric
      stat <- ch.test0(seq_len(S-1))
      pval <- switch(pvalue,
        "raw" = uroot.raw.pvalue(stat, "CH", S-1),
        "RS" = ch.rs.pvalue(stat, type, lag1, S, Nc, rs.nobsreg, S-1))
      stat <- cbind("joint" = stat, "p-value" = pval)
    }
  } else {
    # "sid" is a numeric vector indicating the index of the 
    # seasonal dummy(s) or cycle(s) to be tested
    if (pvalue == "RS")
    {
      #it could be checked if a case for which "RS" p-values are available 
      #is requested through "sid", but for that end, simply use default "sid" value
      pvalue <- "raw"
      warning("argument ", sQuote("pvalue"), " was changed to ", sQuote("raw"))
    }

    #if length(id) > 1, a joint test for those seasons or cycles is obtained   
    stat <- ch.test0(id)
    pval <- uroot.raw.pvalue(stat, "CH", length(id))
    stat <- rbind(c(stat, pval))
    colnames(stat) <- c("statistic", "p-value")   
    if (type == "dummy")
    {
      rownames(stat) <- switch(as.character(S), 
        "4" = paste("Quarter(s)", paste0(id, collapse = ",")), 
        "12" = paste(month.abb[id], collapse = ","),
        paste("Season(s)", paste0(id, collapse = ",")))
    } else { # type trigonometric
##TODO S odd (see how this is arranged above)
    }
  }

  # output

  res <- list(statistics = stat[,1], pvalues = stat[,2], 
    method = "Canova and Hansen test for seasonal stability", data.name = data.name, 
    type = type, fitted.model = fit, 
    NW.order = NW.order, lag1 = lag1, isNullxreg = isNullxreg, type.pvalue = pvalue,
    pvlabels = symnum(stat[,"p-value"], corr = FALSE, na = FALSE,
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
      symbols   =  c("***","**","*","."," ")))
  class(res) <- "CHtest"
  res
}
