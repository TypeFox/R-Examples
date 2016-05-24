#' Order variables in a MLM according to canonical structure or other criteria

varOrder <- function(x, ...) {
	UseMethod("varOrder")
}

varOrder.mlm <- function(x,
                     term,
                     variables, 
                     type = c("can", "pc"),
                     method = c("angles", "dim1", "dim2", "alphabet", "data", "colmean"),
                     names = FALSE,
                     descending = FALSE,
                     ...)
{
	data <- model.frame(x)
	Y <- model.response(data) 
	vars <- colnames(Y)
  if (!missing(variables)){
      if (is.numeric(variables)) {
          vars <- vars[variables]
          if (any(is.na(vars))) stop("Bad response variable selection.")
          }
      else {
          check <- !(variables %in% vars)
          if (any(check)) stop(paste("The following", 
              if (sum(check) > 1) "variables are" else "variable is",
              "not in the model:", paste(variables[check], collapse=", ")))
          vars <- variables
          }
      }

	method = match.arg(method)
	type   = match.arg(type)

	if (method %in% c("angles", "dim1", "dim2")) {
  	if (type == "pc")
  		struc <- eigen(cor(Y))$vectors
  	else
  		struc <- candisc(x, term)$structure
	}
	order <- switch( method,
		alphabet = order(vars),
		angles = order( ifelse( struc[vars,1] >0, 
		                        atan(struc[vars,2]/struc[vars,1]), 
		                        atan(struc[vars,2]/struc[vars,1]) + pi)),
		dim1 = order(struc[vars,1]),
		dim2 = order(struc[vars,2]),
		data = seq_along(vars),
		colmean = order(colMeans(Y))
	)
  
  if (descending) order <- rev(order)
	if (names) vars[order] else order
}

varOrder.data.frame <- 
	function(x, 
           variables, 
           method = c("angles", "dim1", "dim2", "alphabet", "data", "colmean"),
           names = FALSE,
	         descending = FALSE,
	         ...) {
	Y <- x
	vars <- colnames(Y)
  if (!missing(variables)){
      if (is.numeric(variables)) {
          vars <- vars[variables]
          if (any(is.na(vars))) stop("Bad response variable selection.")
          }
      else {
          check <- !(variables %in% vars)
          if (any(check)) stop(paste("The following", 
              if (sum(check) > 1) "variables are" else "variable is",
              "not in the data:", paste(variables[check], collapse=", ")))
          vars <- variables
          }
      }

	method = match.arg(method)
	if (method %in% c("angles", "dim1", "dim2")) {
	  struc <- eigen(cor(Y))$vectors
	}
	order <- switch( method,
  	alphabet = order(vars),
  	angles = order( ifelse( struc[vars,1] >0, 
  	                        atan(struc[vars,2]/struc[vars,1]), 
  	                        atan(struc[vars,2]/struc[vars,1]) + pi)),
  	dim1 = order(struc[vars,1]),
  	dim2 = order(struc[vars,2]),
  	data = seq_along(vars),
  	colmean = order(colMeans(Y))
	)
	if (descending) order <- rev(order)
	if (names) vars[order] else order

}
varOrder.default <- function(x, ...) {
	stop("no methods are yet available for objects of class ", class(x))
}
