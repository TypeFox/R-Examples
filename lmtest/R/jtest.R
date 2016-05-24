## see Greene (2003), Section 8.3.3., p.154

jtest <- function(formula1, formula2, data = list(), vcov. = NULL, ...)
{
  ## merge two models (if possible) and compute
  ## response y and regressor matrices X and Z
  ## 1. If formulas are available: build big model first
  if(inherits(formula1, "formula") && inherits(formula2, "formula")) {
    formula <- formula2
    if(length(formula) > 2) formula[[2]] <- NULL
    formula[[2]] <- as.call(list(as.name("+"), as.name("."), formula[[2]]))
    formula <- update(formula1, formula)
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(terms(formula1), mf)
    Z <- model.matrix(terms(formula2), mf)
    m1 <- deparse(formula1)
    m2 <- deparse(formula2)
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
    if(!identical(y, y2)) warning("different dependent variables specified")
  }
  ## for pretty printing
  m1 <- paste(m1, collapse = "\n")
  m2 <- paste(m2, collapse = "\n")

  ## check vcov.
  if(!is.null(vcov.) && !is.function(vcov.)) stop("`vcov.' needs to be a function")

  ## fit auxiliary models:
  y.hat1 <- lm.fit(X, y)$fitted
  y.hat2 <- lm.fit(Z, y)$fitted
  auxX <- cbind(X, y.hat2)
  auxZ <- cbind(Z, y.hat1)
  auxmodX <- lm(y ~ 0 + auxX)
  auxmodZ <- lm(y ~ 0 + auxZ)

  rvalX <- coeftest(auxmodX, vcov. = vcov., ...)
  rvalZ <- coeftest(auxmodZ, vcov. = vcov., ...)
  rvalX <- rvalX[nrow(rvalX), ]
  rvalZ <- rvalZ[nrow(rvalZ), ]
  
  ## put results together
  title <- "J test\n"
  topnote <- paste("Model ", 1:2,": ", c(m1, m2), sep="", collapse="\n")
  rval <- rbind(rvalX, rvalZ)
  rownames(rval) <- c("M1 + fitted(M2)", "M2 + fitted(M1)")
  rval <- structure(as.data.frame(rval), heading = c(title, topnote),
	            class = c("anova", "data.frame"))
  return(rval)
}
