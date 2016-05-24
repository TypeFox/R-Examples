#' Embed two matrices into one blockdiagonal matrix
#' 
#' @param M matrix of type character
#' @param N matrix of type character
#' @return Matrix of type character containing M and N as upper left and lower right block
#' @export
blockdiagSymb <- function(M, N) {
  
  red <- sapply(list(M, N), is.null)
  if(all(red)) {
    return()
  } else if(red[1]) {
    return(N)
  } else if(red[2]) {
    return(M)
  }
  
  A <- matrix(0, ncol=dim(N)[2], nrow=dim(M)[1])
  B <- matrix(0, ncol=dim(M)[2], nrow=dim(N)[1])
  result <- rbind(cbind(M, A), cbind(B, N))
  return(result)
  
}



#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out data.frame or matrix or list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export
wide2long <- function(out, keep, na.rm) {
  
  UseMethod("wide2long", out)
  
  
}

#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out data.frame or matrix or list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export wide2long.data.frame
#' @export
wide2long.data.frame <- function(out, keep = 1, na.rm = FALSE) {
  
  wide2long.matrix(out, keep = keep, na.rm = na.rm)
  
}

#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out data.frame or matrix or list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export wide2long.matrix
#' @export
wide2long.matrix <- function(out, keep = 1, na.rm = FALSE) {
  
  timenames <- colnames(out)[keep]
  allnames <- colnames(out)[-keep]
  times <- out[,keep]
  ntimes<- nrow(out)
  values <- unlist(out[,allnames])
  outlong <- data.frame(times, name = rep(allnames, each=ntimes), value = as.numeric(values))
  colnames(outlong)[1:length(keep)] <- timenames
  
  if(na.rm) outlong <- outlong[!is.na(outlong$value),]
  
  return(outlong)
  
}

#' Translate wide output format (e.g. from ode) into long format 
#' 
#' @param out list of matrices in wide format 
#' @param keep Index vector, the columns to keep
#' @param na.rm Logical, if \code{TRUE}, missing values are removed in the long format.
#' @details The function assumes that out[,1] represents a time-like vector
#' whereas out[,-1] represents the values. Useful for plotting with ggplot. If 
#' a list is supplied, the names of the list are added as extra column names "condition"
#' @return data.frame in long format, i.e. columns "time" (out[,1]), "name" (colnames(out[,-1])), 
#' "value" (out[,-1]) and, if out was a list, "condition" (names(out))
#' @export wide2long.list
#' @export
wide2long.list <- function(out, keep = 1, na.rm = FALSE) {
  
  conditions <- names(out)
  
  outlong <- do.call(rbind, lapply(conditions, function(cond) {
    
    cbind(wide2long.matrix(out[[cond]]), condition = cond)
    
    #myout <- out[[cond]]
    #timename <- colnames(myout)[1]
    #allnames <- colnames(myout)[-1]
    #times <- myout[,1]
    #values <- unlist(myout[,allnames])
    #myoutlong <- data.frame(time = times, 
    #                        name = rep(allnames, each=length(times)), 
    #                        value = as.numeric(values), 
    #                        condition = cond)
    #colnames(myoutlong)[1] <- timename
    #return(myoutlong)
    
  }))
  
  
  
  return(outlong)
  
}


#' Translate long to wide format (inverse of wide2long.matrix) 
#' 
#' @param out data.frame in long format 
#' @return data.frame in wide format 
#' @export
long2wide <- function(out) {
  
  timename <- colnames(out)[1]
  times <- unique(out[,1])
  allnames <- unique(as.character(out[,2]))
  M <- matrix(out[,3], nrow=length(times), ncol=length(allnames))
  M <- cbind(times, M)
  colnames(M) <- c(timename, allnames)
  
  return(M)
  
}


#' Bind named list of data.frames into one data.frame
#' 
#' @param mylist A named list of data.frame. The data.frames are expected to have the same structure.
#' @details Each data.frame ist augented by a "condition" column containing the name attributed of
#' the list entry. Subsequently, the augmented data.frames are bound together by \code{rbind}.
#' @return data.frame with the originial columns augmented by a "condition" column.
#' @export
lbind <- function(mylist) {
  
  conditions <- names(mylist)
  
  outlong <- do.call(rbind, lapply(conditions, function(cond) {
    
    myout <- mylist[[cond]]
    myoutlong <- cbind(myout, condition = cond)
    
    return(myoutlong)
    
  }))
  
  return(outlong)
  
}

#' Alternative version of expand.grid
#' @param seq1 Vector, numeric or character
#' @param seq2 Vector, numeric or character
#' @return Matrix ob combinations of elemens of \code{seq1} and \code{seq2}
expand.grid.alt <- function(seq1, seq2) {
  cbind(Var1=rep.int(seq1, length(seq2)), Var2=rep(seq2, each=length(seq1)))
}

#' Load a template file in the editor
#' 
#' @param i Integer, choose a template to be loaded
#' @details Possible templates are:
#' i = 1: Do parameter estimation in a dynamic model with fixed forcings
#' @export
loadTemplate <- function(i = 1) {
  
  path <- path.package("dMod")
  if(i == 1) {
    system(paste0("cp ", path, "/templates/R2CTemplate.R mymodel.R"))
    file.edit("mymodel.R")
  }
  
}





#' Evaluation of algebraic expressions defined by characters
#' 
#' @param x Name character vector, the algebraic expressions
#' @param compile Logical. The function is either compiled (requires the \code{inline} package) or
#' evaluated in raw R.
#' @return A prediction function \code{f(mylist)} where \code{mylist} is a list of numeric vectors that can
#' be coerced into a matrix. The names correspond to the symbols used in the algebraic expressions. The
#' function \code{f} returns a matrix.
#' @examples 
#' \dontrun{
#' myfun <- funC0(c(x = "x", y = "a*x^4 + b*x^2 + c"))
#' out <- myfun(list(a = -1, b = 2, c = 3, x = seq(-2, 2, .1)))
#' plot(out[, 1], out[, 2])
#' }
#' 
#' @export
funC0 <- function(x, compile = TRUE) {
    
  # Get symbols to be substituted by x[] and y[]
  outnames <- names(x)
  innames <- getSymbols(x)
  
  x.new <- paste0(x, collapse = ", ")
  x.new <- paste0("list(", x.new, ")")
  x.expr <- parse(text = x.new)
  
  ## Compiled version based on inline package
  ## Non-compiled version based on with() and eval()
  if(compile) {
    
    # Do the replacement to obtain C syntax
    x <- replaceOperation("^", "pow", x)
    x <- replaceSymbols(innames, paste0("x[", (1:length(innames))-1, "+i* *k]"), x)
    names(x) <- paste0("y[", (1:length(outnames)) - 1, "+i* *l]")
    
    # Paste into equation
    x <- x[x != "0"]
    expr <- paste(names(x), "=", x, ";")
    
    # Put equation into loop, body of the C function
    body <- paste(
      "for(int i = 0; i< *n; i++) {",
      paste(expr, collapse=""),
      "}"
    )
    
    # Generate the C function by the inline package
    myCfun <- inline::cfunction(sig=c(x = "double", y = "double", n = "integer", k = "integer", l = "integer"),
                                body=body,
                                language="C",
                                convention=".C"
    )
    
    # Generate output function
    myRfun <- function(x) {
      
      # Translate the list into matrix and then into vector
      M <- do.call(rbind, x[innames])
      if(length(M) == 0) M <- matrix(0)
      x <- as.double(as.vector(M))
      
      # Get integers for the array sizes
      n <- as.integer(dim(M)[2])
      k <- as.integer(length(innames))
      if(length(k) == 0) k <- as.integer(0)
      l <- as.integer(length(outnames))
      
      
      # Initialize output vector
      y <- double(l*n)
      
      # Evaluate C function and write into matrix
      out <- matrix(myCfun(x, y, n, k, l)$y, nrow=length(outnames), ncol=n)
      rownames(out) <- outnames
      
      return(t(out))    
      
    }
    
    
    
  } else {
    
    # Generate output function
    myRfun <- function(x) {
      
      out.list <- with(x, eval(x.expr))
      out.matrix <- do.call(cbind, out.list)
      colnames(out.matrix) <- outnames
      rownames(out.matrix) <- NULL
      
      return(out.matrix)
      
    }
    
  }
  
  
  attr(myRfun, "equations") <- x
  
  return(myRfun)
  
}

