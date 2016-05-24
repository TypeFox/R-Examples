"coca.formula" <-
  function(formula, data, method = c("predictive", "symmetric"),
           reg.method = c("simpls", "eigen"), weights = NULL,
           n.axes = NULL, symmetric = FALSE, ...)
  {
    parseFormula <- function (formula, data)
      {
        Terms <- terms(formula, "Condition", data = data)
        flapart <- fla <- formula <- formula(Terms, width.cutoff = 500)
        specdata <- formula[[2]]
        Yresponse <- as.matrix(eval(specdata, data, parent.frame()))
        formula[[2]] <- NULL
        Ypredictors <- eval(formula[[2]], data, parent.frame())
        #if(class(Ypredictors) == "data.frame")
        #  {
        #    return(list(Yresponse = Yresponse, Ypredictors = Ypredictors))
        #  } else {
            if (formula[[2]] == "1" || formula[[2]] == "0")
              Ypredictors <- NULL
            else {
              mf <- model.frame(formula, data, na.action = na.fail)
              Ypredictors <- model.matrix(formula, mf)
              if (any(colnames(Ypredictors) == "(Intercept)")) {
                xint <- which(colnames(Ypredictors) == "(Intercept)")
                Ypredictors <- Ypredictors[, -xint, drop = FALSE]
              }
            }
        #}
        list(Yresponse = Yresponse, Ypredictors = Ypredictors)
      }
    if (missing(data))
      data <- parent.frame()
    dat <- parseFormula(formula, data)
    x <- dat$Ypredictors
    y <- dat$Yresponse
    nam.dat <- list(namY = deparse(formula[[2]], width.cutoff = 500),
                    namX = deparse(formula[[3]], width.cutoff = 500))
    if (nam.dat$namX == ".")
      nam.dat$namX <- deparse(substitute(data))
    if(any(rowSums(y) <= 0 ))
      stop("all row sums must be >0 in data matrix y")
    if(any((csum <- colSums(y)) <= 0 )) {
      y <- y[, csum > 0, drop = FALSE]
      message("some species contain no data and were removed from data matrix y\n")
    }
    if(any(rowSums(x) <= 0 ))
      stop("all row sums must be >0 in data matrix x")
    if(any((csum <- colSums(x)) <= 0 )) {
      x <- x[, csum > 0, drop = FALSE]
      message("some species contain no data and were removed from data matrix x\n")
    }
    method <- match.arg(method)
    if(method == "predictive")
      {
        reg.method <- match.arg(reg.method)
        retval <- switch(reg.method,
                           simpls = predcoca.simpls(y, x, R0 = weights,
                             n.axes = n.axes, nam.dat),
                           eigen = predcoca.eigen(y, x, R0 = weights,
                             n.axes = n.axes, nam.dat))
      } else {
        retval <- symcoca(y, x, n.axes = n.axes, R0 = weights,
                          symmetric = symmetric, nam.dat)
      }
    retval
  }

