#' Expand an edit expression
#'
#' Often many numeric variables have the same constraints. \code{expandEdits} is
#' a utility function to define edits for multiple variables. See the examples for the syntax.
#' @param s edit expression, can be a \code{character} or \code{expression} vector
#' @param prefix prefix for variables to be expanded
#' @param useSum if \code{TRUE} sum expressions will be expanded
#' @param asExpression if \code{TRUE} an \code{\link{expression}} will be returned in stead of a \code{character}
#' @param env enviroment that will be used to find variables to be expanded
#' @param ... variables used in the expansion
#' @return \code{character} or \code{expression} vector with expanded expressions
#' @keywords internal
expandEdits <- function(s, prefix="_", useSum=TRUE, asExpression=is.language(s), env=parent.frame(), ...){  
  #TODO replace special regex character in prefix with escaped character.
  
  force(asExpression)

  if (is.expression(s)){
    s <- as.character(s)
  }
  
  if (is.language(s)){
    s <- deparse(s)
  }
  
  if (length(s) > 1){
    return(lapply(s, expandEdits, prefix=prefix, useSum=useSum, ...))
  }
  
  l <- as.list(env)
  vars <- list(...)
  l[names(vars)] <- vars
  
  varnms <- paste(prefix,names(l), sep="")
  used <- sapply(varnms, grepl, s)
  varnms <- varnms[used]
  l <- l[used]
  
  if (useSum) {
    sumnms <- paste("\\bsum", names(l), sep=prefix)
    sumregex1 <- paste(sumnms, "\\((.+?)\\)", sep="")
    sumregex2 <- paste(sumnms, "\\((.+?)\\).+", sep="")
    
    vars <- names(l)
    for (i in seq_along(vars)){      
      if (length(grep(sumregex2[i], s))){
        sumvars <- sub(sumregex2[i], "\\1", s)
        sumvars <- do.call(expandEdits, list(s=sumvars, env=l[vars[i]], prefix=prefix))
        sumvars <- paste(sumvars, collapse=" + ")
        s <- sub(sumregex1[i], sumvars, s)
        l[[vars[i]]] <- NULL
      }
    }
  }
  
  varnms <- paste(prefix,names(l), sep="")
  for (i in seq_along(l)){
    if (length(grep(varnms[i], s))) {
      s <- sapply(l[[i]], function(j) gsub(varnms[i],j,s))
    }
  }
  
  if (is.array(s)) {
    dimnames(s) <- l
  } else if (is.vector(s) && length(l)){
    names(s) <- l[[1]]
  }
  if (asExpression){
    parse(text=s)
  } else {
    s
  }
}

