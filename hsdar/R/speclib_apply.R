# if (!isGeneric("apply")) {
#   setGeneric("apply", function(X, MARGIN, FUN, ...)
#   standardGeneric("apply"))
# }

setMethod("apply", signature(X = "Speclib"), 
          function(X, 
                   FUN,
                   byattributes = NULL,
                   ...)
{ 
  usage <- usagehistory(X)
  x <- X 

  if (is.character(FUN))
  {
    FUN_str <- FUN
  } else {
    call_fu <- match.call(call = sys.call(-1))
    FUN_str <- as.character(call_fu[which(names(call_fu) == "FUN")])
  }
  
  if (!is.function(try(match.fun(FUN), silent = TRUE)))
    stop("Unknown function")
  
  FUN <- match.fun(FUN_str)
  
  if (is.null(byattributes))
  {
    result <- speclib(spectra = matrix(data = 0,
                                        ncol = length(wavelength(x)),
                                        nrow = 1),  
                      wavelength = wavelength(x),                      
                      usagehistory = x@usagehistory
                    )
    idSpeclib(result) <- FUN_str  
    
    spectra <- t(spectra(x))
    x <- wavelength(x)
    n <- ncol(spectra)
    
    spec <- apply(spectra, 1, FUN = FUN, ...)
    spec <- matrix(spec, nrow = 1)
    spectra(result) <- spec
    idSpeclib(result) <- FUN_str
    usagehistory(result) <- paste(X@ylabel, "=", FUN_str, "applied to matrix of",n , "spectra")
  } else {
    attributes_X <- attribute(X)
    attributes_col <- which(names(attributes_X) == byattributes)
    if (length(attributes_col) == 0)
      stop(paste("Could not find column '", byattributes, "' in attributes of X", sep = ""))
    attributes_X <- attributes_X[, attributes_col]
    if (!is.factor(attributes_X))
      attributes_X <- as.factor(attributes_X)
    lev <- levels(attributes_X)
    tmp <- parse(text = byattributes)
    tmp_lev <- lev[1]
    result <- apply(.subset.speclib(X, expression(eval(parse(text = tmp)) == tmp_lev)), FUN = FUN, ...)
    attribute(result) <- data.frame(X = lev[1], stringsAsFactors = FALSE)
    for (i in 2:length(lev))
    {
      tmp_lev <- lev[i]
      res_lev <- apply(.subset.speclib(X, expression(eval(parse(text = tmp)) == tmp_lev)),
                       FUN = FUN, ...)
      attribute(res_lev) <- data.frame(X = lev[i], stringsAsFactors = FALSE)
      result <- merge(result, res_lev)
    }
    
    usagehistory(result) <- NULL
    usagehistory(result) <- usage
    usagehistory(result) <- paste(X@ylabel, " = ", FUN_str, " applied to matrix spectra by attribute '", 
                                  byattributes, "'", sep = "")
    names(attribute(result)) <- byattributes
  }
  return(result)
}
)
   




