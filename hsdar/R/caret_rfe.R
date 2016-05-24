if (!isGeneric("rfe")) {
  setGeneric("rfe")
}


setMethod("rfe", signature(x = "Speclib"),
          definition = function(x,
                                y,
                                cutoff = .95,
                                returnData = TRUE,
                                ...)
{  
  y_missing <- missing(y)
  
  if (y_missing)
  {
    y <- .getResponseVar(x, 
                         advice = c("rfe", "setResponse", 
                                    "This is only required if you do not specify 'y'."))
  }
  
  useAttributesAsPredicants <- !is.na(.getPredicantVar(x, stopifmissing = FALSE))[1]
  
  x_dat <- as.data.frame(spectra(x))
  if (is.finite(cutoff))
  {
    x_dat <- x_dat[, -findCorrelation(cor(x_dat), cutoff)]
    x_dat <- as.data.frame(x_dat)
  }
  
  spec_nam <- names(x_dat)
  
  if (useAttributesAsPredicants)
  {
    addVar <- .getPredicantVar(x)  
    x_dat <- cbind(x_dat, addVar)
    if (nlevels(as.factor(names(x_dat))) != ncol(x_dat))
    {
      print(names(x_dat))
      stop("Names in predictor data.frame not unique")
    }
  }

  dots <- list(...)
  res <- if (!any(names(dots) == "rfeControl"))
           rfe(x_dat, y, rfeControl = rfeControl(functions = rfFuncs), ...)
         else
           rfe(x_dat, y, ...)
  if (!returnData)
    return(res)
  
  pred <- predictors(res)
  
  x <- x[,sapply(spec_nam, FUN = function(x, pred) any(pred == x), pred), usagehistory = FALSE]
  
  
  if (useAttributesAsPredicants)
  {
    warning(paste("Attibute data.frame will only contain relevant variables", 
                  if (y_missing) " and the response variable", ".", sep = ""))
    if (y_missing)
      pred <- c(pred, names(attribute(x))[.getCaretParameter(x, "response")])
    cols_keep <- sapply(names(attribute(x)), FUN = function(x, pred) any(pred == x), pred)
    if (sum(cols_keep) > 0)
    {
      if (sum(cols_keep) == 1)
      {
        tmp <- as.data.frame(matrix(attribute(x)[,cols_keep], ncol = 1))
        names(tmp) <- names(attribute(x))[cols_keep]
      } else {
        tmp <- attribute(x)[,sapply(names(attribute(x)), FUN = function(x, pred) any(pred == x), pred)]
      }
      attribute(x) <- tmp
    }
    x <- .updateCaretParameters(x, c("response", "predictor"))
  }
  
  x <- .setCaretParameter(x, "rfe_result", res)
  usagehistory(x) <- "Recursive feature selection"
  return(x)
})

setMethod("rfe", signature(x = "Nri"),
          definition = function(x,
                                y,
                                cutoff = .95,
                                returnData = TRUE,
                                ...)
{ 
  y_missing <- missing(y)
  
  if (y_missing)
  {
    y <- .getResponseVar(x,
                         advice = c("rfe", "setResponse", 
                                    "This is only required if you do not specify 'y'."))
  }  
    
  useAttributesAsPredicants <- !is.na(.getPredicantVar(x, stopifmissing = FALSE))[1]
  
  nri_vals_all <- as.data.frame(x)
  if (is.finite(cutoff))
  {
    nri_vals <- nri_vals_all[, -findCorrelation(cor(nri_vals_all), cutoff)]
  } else {
    nri_vals <- nri_vals_all
  }  
  nri_vals <- as.data.frame(nri_vals)
  
  if (useAttributesAsPredicants)
  {
    addVar <- .getPredicantVar(x)
    nri_vals <- cbind(nri_vals, addVar)
    if (nlevels(as.factor(names(nri_vals))) != ncol(nri_vals))
    {
      print(names(nri_vals))
      stop("Names in predictor data.frame not unique")
    }
  }

  dots <- list(...)
  res <- if (!any(names(dots) == "rfeControl"))
           rfe(nri_vals, y, rfeControl = rfeControl(functions = rfFuncs), ...)
         else
           rfe(nri_vals, y, ...)

  if (!returnData)
    return(res)
  
  pred <- predictors(res)
  
  is.pred.col <- sapply(names(nri_vals_all), FUN = function(x, pred) any(pred == x), pred)
  
  values <- numeric(length = length(x@nri@values))
  values[] <- NA
  incr <- length(x@nri@values)/nrow(nri_vals)
  for (i in 1:ncol(nri_vals_all))
  { 
    if (is.pred.col[i])
    {
      index <- seq(i, length(values), incr)
      values[index] <- nri_vals_all[,i]
    }
  }
  
  x@nri <- distMat3D(values, ncol = ncol(x@nri), nlyr = nrow(nri_vals))  

  if (useAttributesAsPredicants)
  {
    warning(paste("Attibute data.frame will only contain relevant variables", 
                  if (y_missing) " and the response variable", ".", sep = ""))
    if (y_missing)
      pred <- c(pred, names(attribute(x))[.getCaretParameter(x, "response")])
    cols_keep <- sapply(names(attribute(x)), FUN = function(x, pred) any(pred == x), pred)
    if (sum(cols_keep) > 0)
    {
      if (sum(cols_keep) == 1)
      {
        tmp <- as.data.frame(matrix(attribute(x)[,cols_keep], ncol = 1))
        names(tmp) <- names(attribute(x))[cols_keep]
      } else {
        tmp <- attribute(x)[,sapply(names(attribute(x)), FUN = function(x, pred) any(pred == x), pred)]
      }
      attribute(x) <- tmp
    }
    x <- .updateCaretParameters(x, c("response", "predictor"))
  }
  
  return(.setCaretParameter(x, "rfe_result", res))
})

setMethod("rfe", signature(x = "Specfeat"),
          definition = function(x,
                                y,
                                cutoff = .95,
                                returnData = TRUE,
                                ...)
{
  x <- .as.speclib.specfeat(x)
  if (missing(y))
  {
    return(rfe(x, cutoff = cutoff, returnData = returnData, ...))
  } else {
    return(rfe(x, y, cutoff = cutoff, returnData = returnData, ...))
  }
})

get_rfe  <- function(x)
  .getCaretParameter(x, "rfe_result")