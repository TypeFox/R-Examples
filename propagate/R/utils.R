########## visible ########################
makeGrad <- function(expr, order = NULL)
{
  VARS <- all.vars(expr)
  if (!is.null(order)) VARS <- VARS[order]
  FUN <- function(x) D(expr, x)
  vecGRAD <- sapply(VARS, FUN)
  vecGRAD <- matrix(vecGRAD, nrow = 1)    
  return(vecGRAD)  
} 

makeHess <- function(expr, order = NULL)
{
  VARS <- all.vars(expr)  
  if (!is.null(order)) VARS <- VARS[order]
  GRID <- expand.grid(VARS, VARS)    
  FUN <- function(x) D(D(expr, x[1]), x[2])
  vecHESS <- apply(GRID, 1, FUN)  
  matHESS <- matrix(vecHESS, ncol = length(VARS), byrow = TRUE)    
  return(matHESS)
} 

evalDerivs <- function(deriv, envir)
{
  if (missing(envir)) envir <- .GlobalEnv
  DIM <- dim(deriv)
  evalVEC <- sapply(deriv, eval, envir = envir)
  dim(evalVEC) <- DIM
  return(evalVEC)
}

kurtosis <- function (x, na.rm = FALSE) 
{
  if (is.matrix(x)) 
    apply(x, 2, kurtosis, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    n <- length(x)
    n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2) - 3
  }
  else if (is.data.frame(x)) 
    sapply(x, kurtosis, na.rm = na.rm)
  else kurtosis(as.vector(x), na.rm = na.rm)
}

skewness <- function (x, na.rm = FALSE) 
{
  if (is.matrix(x)) 
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x)) 
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

counter <- function (i) 
{
  if (i%%10 == 0) 
    cat(i)
  else cat(".")
  if (i%%50 == 0) 
    cat("\n")
  flush.console()
}

tr <- function(mat) sum(diag(mat), na.rm = TRUE)

rescale <- function (x, tomin, tomax) 
{
  if (missing(x) | missing(tomin) | missing(tomax)) {
    stop(paste("Usage: rescale(x, tomin, tomax)\n", "\twhere x is a numeric object and tomin and tomax\n is the range to rescale into", 
               sep = "", collapse = ""))
  }
  if (is.numeric(x) && is.numeric(tomin) && is.numeric(tomax)) {
    xrange <- range(x, na.rm = TRUE)
    if (xrange[1] == xrange[2]) 
      return(x)
    mfac <- (tomax - tomin)/(xrange[2] - xrange[1])
    return(tomin + (x - xrange[1]) * mfac)
  }
  else {
    warning("Only numeric objects can be rescaled")
    return(x)
  }
}

print.interval <- function(x, ...)
{
  cat("[", x[1], ", ", x[2], "]", sep = "")
}

