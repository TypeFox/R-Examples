kmData <- function(formula, data, inputnames=NULL, ...) {
	  
  f <- as.formula(formula)
  data <- as.data.frame(data)  # in case where data is a matrix
  tt <- terms(x=f, data = data)
  
  if (attr(tt, "response") == 0) {
    stop("the response must be provided in formula")
  } else {
    yname <- all.vars(f, max.names = 1)
    datanames <- names(data)
    if (!(yname %in% datanames)) { 
      stop("the response name (left hand side of the formula) must belong to the variable names specified in 'data'")
    } else {
      ypos <- match(yname, datanames)
      response <- data[, ypos, drop=FALSE]
      xnames <- setdiff(datanames, yname)
      if (is.null(inputnames)) {
        inputnames <- xnames
      }
      if (length(intersect(inputnames, xnames))!=length(inputnames)) {
        stop("the input names must belong to the variable names specified in 'data' (and do not contain the response name)")  
      } else {
        xpos <- match(inputnames, datanames)
        design <- data[, xpos, drop=FALSE]
      }      
    }  
    
    m <- km(formula=formula, response=response, design=design, ...)
    return(m)
  }

}
