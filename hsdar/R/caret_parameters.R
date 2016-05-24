if (!isGeneric("setResponse")) {
  setGeneric("setResponse", function(x, response)
  standardGeneric("setResponse"))
}

if (!isGeneric("setPredictor")) {
  setGeneric("setPredictor", function(x, predictor)
  standardGeneric("setPredictor"))
}

if (!isGeneric("showCaretParameters")) {
  setGeneric("showCaretParameters", function(x)
  standardGeneric("showCaretParameters"))
}

setMethod("setResponse", signature(x = ".CaretHyperspectral", response = "character"),
          definition = function(x, response)
{
  response_index <- sapply(response, function(response, x)
    {
      ind <- which(x == response)
      if (length(ind) == 0)
        stop(paste("'", response, "' not found in ", class(x), " x", sep = ""))
      return(ind)
    }, names(attribute(x)))
  
  x <- .setCaretParameter(x, "response", response_index)
  x <- .setCaretParameter(x, "responseName", response)
  usagehistory(x) <- paste0("Response variable(s) set to \"", paste0(response, collapse = "\", \""), "\"")
  return(x)
})


setMethod("setPredictor", signature(x = ".CaretHyperspectral", predictor = "character"),
          definition = function(x, predictor)
{
  
  predictor_index <- sapply(predictor, function(predictor, x)
    {
      ind <- which(x == predictor)
      if (length(ind) == 0)
        stop(paste("'", predictor, "' not found in ", class(x), " x", sep = ""))
      return(ind)
    }, names(attribute(x)))
  
  x <- .setCaretParameter(x, "predictor", predictor_index)
  x <- .setCaretParameter(x, "predictorName", predictor)
  usagehistory(x) <- paste0("Predictor variable(s) set to \"", paste0(predictor, collapse = "\", \""), "\"")
  return(x)
})


.setCaretParameter <- function(x, parameter, value, usagehistory = NULL)
{
  if (is.null(attr(x, "caretParameters"))) ## create new
  {    
    tmp <- list(parameter = value)
    names(tmp) <- parameter
  } else {
    tmp <- attr(x, "caretParameters")
    if (parameter %in% names(tmp)) ## update 
    {
      if (length(value) == 0)
      {
        tmp[[which(parameter == names(tmp))]] <- NA
      } else {
        tmp[[which(parameter == names(tmp))]] <- value
      }
    } else {   ## add
      tmp$parameter <- value
      names(tmp)[length(tmp)] <- parameter
    }    
  }
  attr(x, "caretParameters") <- tmp
  if (!is.null(usagehistory))
  {
    if (is.speclib(x))
      usagehistory(x) <- usagehistory
  }
  return(x)
}

.updateCaretParameters <- function(x, parameters)
{
  
  for (i in 1:length(parameters))
  {
    para <- .getCaretParameter(x, parameters[i])
    paraName <- .getCaretParameter(x, paste(parameters[i], "Name", sep = ""),
                                   stopifmissing = FALSE)
    if (!is.na(paraName[1]))
    {
      if (parameters[i] %in% .getAttrParameters())
      {
        still_valid <- sapply(paraName, function(x, avl)
          {
            ind <- which(avl == x)
            if (length(ind) == 0)
            {
              return(0)
            } else {
              return(ind)
            }
          }, names(attribute(x)))
        x <- eval(parse(text = paste("set", toupper(substr(parameters[i], 1, 1)),
                                     substr(parameters[i], 2, nchar(parameters[i])), 
                                     "(x, paraName[still_valid > 0])", 
                                     sep = "")))
      }
    }
  }
  return(x)
}

.getCaretParameter <- function(x, parameter, advice = NULL, stopifmissing = TRUE)
{
  if (is.null(attr(x, "caretParameters")))
  { 
    if (stopifmissing)
    {
      stop(paste("Object does not contain caretParameters.", 
                 if (!is.null(advice)) paste("Please run function '",advice[2], "' prior to '",
                                             advice[1], "'.", sep = ""),
                 if (!is.null(advice) & length(advice) > 2) paste("\n  ", advice[3], sep = "")))
    } else {
      return(NA)
    }
  }
  tmp <- attr(x, "caretParameters")
  
  if (!(parameter %in% names(tmp)))
  {
    if (stopifmissing)
    {
      stop(paste("Object does not contain required parameter(s)", 
                if (!is.null(advice)) paste("Please run function '",advice[2], "' prior to '",
                                            advice[1], "'.", sep = ""),
                if (!is.null(advice) & length(advice) > 2) paste("\n  ", advice[3], sep = "")))
    } else {
      return(NA)
    }
  }  
  return(tmp[[which(parameter == names(tmp))]])
}

.getPredicantVar <- function(x, ...)
{
  varIndex <- .getCaretParameter(x, "predictor", ...)
  if (is.na(varIndex[1]))
    return(NA)
  attribute(x)[,varIndex]
}

.getResponseVar <- function(x, ...)
{
  varIndex <- .getCaretParameter(x, "response", ...)
  if (is.na(varIndex[1]))
    return(NA)
  attribute(x)[,varIndex]
}

setMethod("showCaretParameters", signature(x = ".CaretHyperspectral"),
          definition = function(x)
{
  if (is.null(attr(x, "caretParameters")))
  {
    cat("Object does not contain caretParameters.\n")
  } else {
    para <- attr(x, "caretParameters")
    for (i in 1:length(para))
    {
      cat(paste("********************************\n",
                names(para)[i], "\n", sep = ""))
      print(para[[i]])
      cat("\n")
    }
  }
})

.getAllPredictors <- function(x, cutoff)
{
  useAttributesAsPredicants <- !is.na(.getPredicantVar(x, stopifmissing = FALSE))[1]

  all_spectra_vals <- as.data.frame(x)
  if (is.finite(cutoff))
  {
    all_vals <- all_spectra_vals[, -findCorrelation(cor(all_spectra_vals), cutoff)]
  } else {
    all_vals <- all_spectra_vals
  }  
  all_vals <- as.data.frame(all_vals)
  
  spectral <- ncol(all_vals)
  if (useAttributesAsPredicants)
  {
    addVar <- .getPredicantVar(x)
    all_vals <- cbind(all_vals, addVar)
  }
  
  attr(all_vals, "spectral") <- c(1:spectral)
  attr(all_vals, "useattributes") <- ncol(all_vals) > spectral
  return(all_vals)
}
