#' An \code{editmatrix} is a numerical matrix and a set of comparison operators representing
#' a linear system of (in)equations.  
#'
#' The function \code{editmatrix} generates an editmatrix from a \code{character} vector, an \code{expression}
#' vector or a \code{data.frame} with at least the column \code{edit}. The function \code{\link{editfile}} 
#' reads edits from a free-form textfile, function \code{\link{as.editmatrix}} converts a matrix, a vector of
#' constants and a vector of operators to an \code{editmatrix}
#'
#' By default, the \code{editmatrix} is normalized, meaning that all comparison operators are converted
#' to one of \code{<}, \code{<=}, or \code{==}. Users may specify edits using any of the operators
#' \code{<, <=, ==, >=, >} (see examples below). However it is highly recommended to let \code{editmatrix} parse them into normal
#' form as all functions operating on editmatrices expect or convert it to normal form anyway.
#'
#' @title Create an editmatrix
#' @seealso 
#'      \code{\link{editrules.plotting}}, \code{\link{violatedEdits}}, \code{\link{localizeErrors}},
#'      \code{\link{normalize}}, \code{\link{contains}}, \code{\link{is.editmatrix}},
#'      \code{\link{getA}}, \code{\link{getAb}}, \code{\link{getb}}, \code{\link{getOps}} \code{\link{getVars}},
#'      \code{\link{eliminate}}, \code{\link{substValue}}, \code{\link{isFeasible}}
#' @export
#' @example ../examples/editmatrix.R
#'
#' @param editrules  A \code{character} or \code{expression} vecotr with (in)equalities written in R syntax.
#'      Alternatively, a \code{data.frame} with a column named \code{edits}, see details.
#' @param normalize \code{logical} specifying if all edits should be transformed (see description)
#'
#' @return \code{editmatrix} : An object of class \code{editmatrix} 
editmatrix <- function( editrules
                      , normalize = TRUE
					       ){   
   if (is.character(editrules) || is.expression(editrules)){
      edit <- editrules
      name <- names(editrules)
      description <- NULL
      editrules <- NULL
    } else if (is.data.frame(editrules)){
      if (is.null(editrules$edit)){
         stop("The supplied data.frame misses the column 'edit'.\nSee ?editmatrix for a valid input specification")
      }            
      name <- editrules$name
      edit <- as.character(editrules$edit)
      description <- editrules$description
    } else {
      stop("Invalid input, please use a character vector, expression vector, or a data.frame.\n See ?editmatrix for a valid input specification")
    }
    if ( length(edit) == 0 ) return(neweditmatrix(matrix(numeric(0)),ops=character(0),normalized=TRUE))
    edts <- parseEdits(edit, type="num")   
  	if (is.null(name)){
  	   name <- paste("num", seq_along(edts),sep="")
  	}
    rowedts <- lapply(edts, function(edt){parseNum(edt)})
    ops <- sapply(edts, function(e){deparse(e[[1]])})
   
    vars <- unique(names(unlist(rowedts)))
    vars <- c(vars[vars!="CONSTANT"], "CONSTANT")

    A <- matrix( 0
                 , ncol=length(vars)
                 , nrow=length(rowedts)
                 , dimnames = list( rules = name
                                  , var=vars
                                  )
                 )
                 
    for (i in 1:length(rowedts)){
       A[i,names(rowedts[[i]])] <- rowedts[[i]]
    }
    A[,ncol(A)] <- -A[,ncol(A)]
   
   if (normalize){
      geq <- ops == ">="
      gt <- ops == ">"
      A[geq | gt,] <- -A[geq | gt,]
      ops[geq] <- "<="
      ops[gt] <- "<"      
   }
   
   names(ops) <- name
   E <- neweditmatrix(A, ops=ops, normalized=all(ops %in% c("==","<","<=")))
   attr(E, "description") <- description
   E
}

#' Create an \code{editmatrix} object from its constituing attributes. 
#'
#' This function is for internal purposes, please use \code{\link{editmatrix}} for creating an editmatrix object.
#' @param A An augmented \code{matrix} of the form \code{A|b}
#' @param ops a character vector with the comparison operator of every edit.
#' @param normalized \code{logical} TRUE or FALSE
#' @param ... optional attributes
#' @return an S3 object of class \code{editmatrix} 
#' @keywords internal
neweditmatrix <- function(A, ops, normalized=FALSE,...){
   structure( A
            , class="editmatrix"
            , ops = ops
            , normalized = normalized
            , ...
            )
}





#' Coerce a matrix to an edit matrix.
#'
#' \code{as.editmatrix} interpretes the matrix as an editmatrix.
#' The columns of the matrix
#' are the variables and the rows are the edit rules (contraints).
#' 
#' If only argument \code{x} is given (the default), the resulting editmatrix is of the form \eqn{Ax=0}. 
#' This can be influenced by using the parameters \code{b} and \code{ops}.
#'
#' @export
#' @seealso \code{\link{editmatrix}}
#'
#' @param A matrix to be transformed into an \code{\link{editmatrix}}. 
#' @param b Constant, a \code{numeric} of \code{length(nrow(x))}, defaults to 0
#' @param ops Operators, \code{character} of \code{length(nrow(x))} with the equality operators, defaults to "=="
#' @param ... further attributes that will be attached to the resulting editmatrix
#'
#' @return an object of class \code{editmatrix}.
as.editmatrix <- function( A
                         , b = numeric(nrow(A))
                         , ops = rep("==", nrow(A))
                         , ...
                         ){
    if (is.editmatrix(A)){
        return(A)
    } 
    if (!is.matrix(A)){
        stop("Argument A must be an object of class matrix. ")
    }   
    cn <- colnames(A)
    if (is.null(cn)){
       cn <- make.names(paste("x", 1:ncol(A), sep=""), unique=TRUE)
    }
    rn <- rownames(A)
    if ( is.null(rn) ){
      rn <- paste("num", 1:nrow(A), sep="")
    } else {
      rn <- make.unique(rn)
    }
    A <- cbind(as.matrix(A), b)
    dimnames(A) <- list(rules=rn,var=c(cn,"CONSTANT"))
    E <- neweditmatrix(A=A, ops=ops, ...)
    if (isNormalized(E)) attr(E,"normalized") <- TRUE
    E
}

#'
#' @note since version 2.0-0, the behaviour of \code{as.data.frame.editmatrix} changed to be more symmetrical 
#' with \code{editmatrix.data.frame} and \code{as.data.frame.editarray}. 
#' Use \code{editrules:::toDataFrame} (unsupported) for the old behaviour.
#'
#' @export 
#' @rdname editmatrix
#' @method as.data.frame editmatrix
#' @param x editmatrix object
#'
#' @return \code{as.data.frame} a 3-column \code{data.frame} with columns 'name' and 'edit'. 
#' If the input editmatrix has a \code{description} attribute a third column is returned.
as.data.frame.editmatrix <- function(x,...){
    edts <- as.character(x,...)
    d <- data.frame(
        name=names(edts),
        edit=edts,
        row.names=NULL,
        stringsAsFactors=FALSE
    )
    if (!is.null(attr(x,'description'))) d$description <- attr(x,'description')
    d
}


# Coerce an editmatrix to a \code{data.frame}
#
# Coerces an editmatrix to a \code{data.frame}. Useful for viewing the matrix representation of editrules.
# 
# 
# @param x \code{\link{editmatrix}} object
#
# @return data.frame with the coefficient matrix representation of \code{x}, an operator column and CONSTANT column.
toDataFrame <- function(x){
   if (!is.editmatrix(x)) stop('x must be an editmatrix')
   dat <- as.data.frame(getA(x))
   nms <- make.names( c(names(dat), "Ops", "CONSTANT")
                    , unique=TRUE
                    )
   n <- length(nms)
   dat[[nms[n-1]]] <- getOps(x)
   dat[[nms[n]]] <- getb(x)
   dat
}



#'
#' @rdname editmatrix
#' @export
#' @method as.character editmatrix
#'
as.character.editmatrix <- function(x, ...){
   A <- getA(x)
   b <- getb(x)
   vars <- getVars(x)
   ops <- getOps(x)

   n <- ncol(A)

   nC <- ncol(A) + 1
   er <- character(nrow(A))

   left <- right <- character(nrow(A)) 
   for ( i in seq_along(rownames(A)) ){
     r <- A[i,]
     lhs <- r > 0
     rhs <- r < 0
     left[i] <- if(any(lhs)) { paste(r[lhs], "*", vars[lhs],sep="",collapse=" + ") } else ""
     right[i] <-if(any(rhs))  { paste(-r[rhs], "*",vars[rhs],sep="",collapse=" + ") } else ""
   }
   left <- gsub(" 1\\*"," ",left)
   right <- gsub(" 1\\*"," ",right)
   left <- gsub("^1\\*","",left)
   right <- gsub("^1\\*","",right)

   
   left <- ifelse( left==""
                 , ifelse(right=="", "0", -b) 
                 , ifelse(b < 0 & right != "", paste(left,-b,sep=" + "), left)
                 )
   
   right <- ifelse( right==""
                  , b
                  , ifelse(b > 0 & left != -b, paste(right,b,sep=" + "), right)
                  )
   txt <- paste(left,ops,right)    
   names(txt) <- rownames(x)
   txt
}

#'
#' @export
#' @rdname editmatrix
#' @method as.expression editmatrix
#'
as.expression.editmatrix <- function(x, ...){
  return(
     tryCatch(parse(text=as.character(x, ...)), 
         error=function(e){
             stop(paste("Not all edits can be parsed, parser returned", e$message,sep="\n"))
         }
     )
 )
}

