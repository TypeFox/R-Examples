NUMCMP <- c("==","<","<=",">",">=")
NUMOPS <- c("+","-","*","<","<=",">=",">")

#' Parse a numerical edit expression 
#'
#' Parse a numerical edit expression into a named \code{numeric}.
#' The \code{names} are the names of the variables
#' @param e a valid R expression
#' @keywords internal
parseNum <- function(e){
  if (!isNum(e)){
     stop(paste("Invalid edit rule:", e))
  }
  wgt <- retrieveCoef(e)
  # simplify the coefficients by summing them
  tapply(wgt, names(wgt), sum)
}


hasNum <- function(e){
  if (length(e) == 1){
    return(is.numeric(e))
  }
  op <- deparse(e[[1]])
  if (length(e) == 2){
    return (op %in% NUMOPS || hasNum(e[[2]]))
  }
  if (length(e) == 3){
    if (op == "%in%") return(FALSE)
    return(op %in% NUMOPS || hasNum(e[[2]]) || hasNum(e[[3]]))
  }
}

# basic test for numerical edit
isNum <- function(e){
  if (length(e) != 3) 
    return(FALSE)  
  cmp <- deparse(e[[1]])
  if (cmp == "==") { 
    return(!is.character(e[[3]]) && !is.logical(e[[3]]))
  }
  return(cmp %in% NUMCMP)
}
  
retrieveCoef <- function(e, co=1){
   #stopifnot(is.language(e))
   if (length(e) == 1){
     if (is.numeric(e)){
        l <- co*e   #the resulting matrix is augmented so b has a -
        names(l) <- getOption("editrules.CONSTANT", "CONSTANT")
     }
     else {
        l <- co
        names(l) <- as.character(e)
     }
     return(l)
   }
   if (length(e) == 2){
     op <- deparse(e[[1]])
      rhs <- e[[2]]
     if (op == "("){
      return(retrieveCoef(rhs, co))
	  } else if (op == "-"){
        return(retrieveCoef(rhs, -1*co))
     }
	  else { 
		stop("Operator ", op, " not implemented", "Invalid expression:", e)
	  }
   }
   if (length(e) == 3){
      op <- deparse(e[[1]])
      lhs <- e[[2]]
      rhs <- e[[3]]
      lsign <- rsign <- co
      if ( op %in% c(NUMCMP, "-")){
	      rsign <- -1 * co
	    } 
	    else if (op == "+"){
	    }
	    else if (op == "*"){
        if (length(lhs) == 2 || is.numeric(lhs)){
          co <- retrieveCoef(lhs, co)
          return(retrieveCoef(rhs, co))          
        } else if (length(rhs) == 2 || is.numeric(rhs)){
          co <- retrieveCoef(rhs, co)
          return(retrieveCoef(lhs, co))         
        } else{
          stop("Expression '", deparse(e), "' contains nonconstant coefficient")
        }
	  }
	  else { 
		   stop("Operator ", op, " not implemented", "Invalid expression:", e)
	  }
	  return(c( retrieveCoef(lhs, lsign)
		        , retrieveCoef(rhs, rsign)
		        )
		      )
   }
   stop("Invalid expression:", e)
}
