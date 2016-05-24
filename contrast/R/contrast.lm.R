## This contrast method is very much like contrast.rms, but uses
## the predictFrame function (defined earlier in this file) instead of
## the predictDesign function defined in the Design package.  It also
## uses the testStatistic function (defined earlier in this file) to make
## it more modular.

## contrast <- function (fit, ...) UseMethod("contrast")

## define the gls, lme and geese versions to execute the lm version
contrast.gls   <- function(fit, ...)
{
  library(nlme)
  contrastCalc(fit, ...)
}

contrast.lme   <- function(fit, ...)
{
  library(nlme)
  contrastCalc(fit, ...)
}

contrast.geese <- function(fit, ...)
{
  library(geepack)
  contrastCalc(fit, ...)
}

contrast.lm    <- function(fit, ...) contrastCalc(fit, ...)


contrastCalc <- function(fit, a, b, cnames=NULL,
                         type=c('individual', 'average'),
                         weights='equal', conf.int=0.95, 
                         fcType = "simple",
                         fcFunc = I,
                         covType = NULL,
                         ..., 
                         env=parent.frame(2))
{
  type <- match.arg(type)

  da <- do.call('generateData', list(fit=fit, factors=a, env=env))
  xa <- predictFrame(fit, da, env=env)
  ma <- nrow(xa)

  if (missing(b))
    {
      xb <- 0 * xa
      db <- da
    } else {
      db <- do.call('generateData', list(fit, factors=b, env=env))
      xb <- predictFrame(fit, db, env=env)
    }
  mb <- nrow(xb)

  vary <- NULL
  if (type == 'individual' && length(cnames) == 0)
    {
      ## If two lists have same length, label contrasts by any variable
      ## that has the same length and values in both lists
      if (ma == mb)
        {
          if (ncol(da) != ncol(db)) stop('program logic error')
          if (any(sort(names(da)) != sort(names(db)))) stop('program logic error')
          k <- integer(0)
          nam <- names(da)
          for (j in 1:length(da))
            {
              if (all(as.character(da[[nam[j]]]) == as.character(db[[nam[j]]])))
                k <- c(k, j)
            }
          if (length(k) > 0) vary <- da[k]
        } else if (max(ma, mb) > 1) {
          ## Label contrasts by values of longest variable in list if
          ## it has the same length as the expanded design matrix
          d <- if (ma > 1) a else b
          l <- sapply(d, length)
          vary <- if (sum(l == max(ma, mb)) == 1) d[l == max(ma, mb)]
        }
    }

  if (max(ma, mb) > 1 && min(ma, mb) == 1)
    {
      if (ma == 1)
        xa <- matrix(xa, nrow=mb, ncol=ncol(xb), byrow=TRUE)
      else
        xb <- matrix(xb, nrow=ma, ncol=ncol(xa), byrow=TRUE)
    } else if (mb != ma) {
      stop('number of rows must be the same for observations generated\nby a and b unless one has one observation')
    }

  X <- xa - xb
  p <- ncol(X)
  m <- nrow(X)

  modelCoef <- getCoefficients(fit)
  denom <- xb %*% modelCoef 
  numer <- xa %*% modelCoef 
  ratio <- fcFunc(numer) / fcFunc(denom) 
  fc <- switch(
               fcType,
               simple = ratio,
               log = log(ratio),
               signed = ifelse(ratio > 1, ratio, -1/ratio))
  
  if (is.character(weights))
    {
      if (weights!='equal') stop('weights must be "equal" or a numeric vector')
      weights <- rep(1, m)
    } else if (length(weights) > 1 && type == 'individual') {
      stop('can specify more than one weight only for type="average"')
    } else if (length(weights) != m) {
      stop(paste('there must be', m, 'weights'))
    }
  weights <- as.vector(weights)

  if (m > 1 && type == 'average')
    X <- matrix(apply(weights * X, 2, sum) / sum(weights), nrow=1,
                dimnames=list(NULL, dimnames(X)[[2]]))

  if(class(fit)[1] == "lm")
    {
      library(sandwich)
      if(is.null(covType)) covType <- "const"
      covMat <- try(vcovHC(fit, type = covType), silent = TRUE)
      if(class(covMat)[1] == "try-error") {
        warning("Sandwich estimate failed; using standard estimate instead")
        covMat <- vcov(fit)
      }
    } else covMat <- vcov(fit)
  
  
  res <- testStatistic(fit, X, modelCoef, covMat, conf.int = conf.int)
  res$cnames <- if (type == 'average') NULL else cnames
  res$nvary <- length(vary)
  res$foldChange <- fc
  res$aCoef <- xa
  res$bCoef <- xb
  res$model <- class(fit)[1]
  res$covType <- covType
  if (type == 'individual') res <- c(vary, res)
  structure(res, class='contrast')
}
