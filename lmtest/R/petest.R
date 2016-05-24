## see Verbeek (2008), Section 3.2.3, p.64

petest <- function(formula1, formula2, data = list(), vcov. = NULL, ...)
{
  ## merge two models (if possible) and compute
  ## response y/log(y) and regressor matrices X and Z
  ## 1. If formulas are available: build big model first
  if(inherits(formula1, "formula") && inherits(formula2, "formula")) {
    ychar1 <- deparse(formula1[[2]])
    ychar2 <- deparse(formula2[[2]])
    if(nchar(ychar1) > nchar(ychar2)) {
      if(substr(ychar1, 5, nchar(ychar1) - 1) != ychar2)
        stop("formula1 does not seem to specify the transformed response of formula2")
      trafo <- substr(ychar1, 1, 3)
      if(!(trafo %in% c("log", "exp")))
        stop("formula1 does not seem to specify the log/exp transformed response of formula2")

      if(trafo == "log") {
        islog1 <- TRUE
        islog2 <- FALSE
      } else {
        islog1 <- FALSE
        islog2 <- TRUE
      }

      formula <- formula1
      if(length(formula) > 2) formula[[2]] <- NULL
      formula[[2]] <- as.call(list(as.name("+"), as.name("."), formula[[2]]))
      formula <- update(formula2, formula)
    } else {
      if(substr(ychar2, 5, nchar(ychar2) - 1) != ychar1)
        stop("formula2 does not seem to specify the transformed response of formula1")
      trafo <- substr(ychar2, 1, 3)
      if(!(substr(ychar2, 1, 3) %in% c("log", "exp")))
        stop("formula2 does not seem to specify the log/exp transformed response of formula1")

      if(trafo == "log") {
        islog1 <- FALSE
        islog2 <- TRUE
      } else {
        islog1 <- TRUE
        islog2 <- FALSE
      }

      formula <- formula2
      if(length(formula) > 2) formula[[2]] <- NULL
      formula[[2]] <- as.call(list(as.name("+"), as.name("."), formula[[2]]))
      formula <- update(formula1, formula)
    }

    mf <- model.frame(formula, data = data)
    yraw <- model.response(mf)
    X <- model.matrix(delete.response(terms(formula1)), mf)
    Z <- model.matrix(delete.response(terms(formula2)), mf)
    m1 <- deparse(formula1)
    m2 <- deparse(formula2)

    if(trafo == "log") {
      y  <- if(islog1) log(yraw) else yraw
      y2 <- if(islog2) log(yraw) else yraw
    } else {
      y  <- if(islog1) yraw else exp(yraw)
      y2 <- if(islog2) yraw else exp(yraw)
    }
  } else {
  ## 2. if not, go via row names
    if(!inherits(formula1, "formula")) {
      mf <- model.frame(formula1)
      X <- if(is.matrix(formula1$x)) formula1$x else model.matrix(terms(formula1), mf)
      y <- if(is.vector(formula1$y)) formula1$y else model.response(mf)
      m1 <- deparse(formula(formula1))
    } else {
      mf <- model.frame(formula1, data = data)
      y <- model.response(mf)
      X <- model.matrix(formula1, data = data)
      m1 <- deparse(formula1)
    }  

    if(!inherits(formula2, "formula")) {
      mf <- model.frame(formula2)
      Z <- if(is.matrix(formula2$x)) formula2$x else model.matrix(terms(formula2), mf)
      y2 <- if(is.vector(formula2$y)) formula2$y else model.response(mf)
      m2 <- deparse(formula(formula2))
    } else {
      mf <- model.frame(formula2, data = data)
      y2 <- model.response(mf)
      Z <- model.matrix(formula2, data = data)
      m2 <- deparse(formula2)
    }  

    if(!(all(c(row.names(X) %in% row.names(Z), row.names(Z) %in% row.names(X))))) {
      warning("models fitted on different subsets")
      allnames <- row.names(X)[row.names(X) %in% row.names(Z)]
      X <- X[allnames,]
      Z <- Z[allnames,]
      y <- y[allnames]
      y2 <- y2[allnames]
    }
    
    if(isTRUE(all.equal(y2, log(y)))) {
      islog1 <- FALSE
      islog2 <- TRUE    
    } else if(isTRUE(all.equal(log(y2), y))) {
      islog1 <- TRUE
      islog2 <- FALSE
    } else {
        stop("response in formula2 does not seem to specify the log or exp of the response in formula1")    
    }
  }
  ## for pretty printing
  m1 <- paste(m1, collapse = "\n")
  m2 <- paste(m2, collapse = "\n")

  ## check vcov.
  if(!is.null(vcov.) && !is.function(vcov.)) stop("`vcov.' needs to be a function")

  ## fit auxiliary models:
  y.hat1 <- lm.fit(X, y)$fitted
  y.hat2 <- lm.fit(Z, y2)$fitted
  y.aux1 <- if(islog1) y.hat2 - exp(y.hat1) else log(y.hat1) - y.hat2
  y.aux2 <- if(islog2) y.hat1 - exp(y.hat2) else log(y.hat2) - y.hat1
  aux1 <- cbind(X, y.aux1)
  aux2 <- cbind(Z, y.aux2)
  auxmod1 <- lm(y ~ 0 + aux1)
  auxmod2 <- lm(y2 ~ 0 + aux2)

  rval1 <- coeftest(auxmod1, vcov. = vcov., ...)
  rval2 <- coeftest(auxmod2, vcov. = vcov., ...)
  rval1 <- rval1[nrow(rval1), ]
  rval2 <- rval2[nrow(rval2), ]
  str1 <- paste("M1 +", if(islog1) "fit(M2)-exp(fit(M1))" else "log(fit(M1))-fit(M2)")
  str2 <- paste("M2 +", if(islog2) "fit(M1)-exp(fit(M2))" else "log(fit(M2))-fit(M1)")
  
  ## put results together
  title <- "PE test\n"
  topnote <- paste("Model ", 1:2,": ", c(m1, m2), sep="", collapse="\n")
  rval <- rbind(rval1, rval2)
  rownames(rval) <- c(str1, str2)
  rval <- structure(as.data.frame(rval), heading = c(title, topnote),
	            class = c("anova", "data.frame"))
  return(rval)
}
