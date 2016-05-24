#' Derive editmatrix with soft constraints based on boundaries of variables. This is a utility function that is used for 
#' constructing a mip/lp problem.
#' @param E normalized \code{editmatrix}
#' @param prefix \code{character} used for naming dummy variables in matrix.
#' @keywords internal
softEdits <- function(E, prefix="delta.", ...){
  UseMethod("softEdits")
}

#' Derive editmatrix with soft constraints based on boundaries of variables. This is a utility function that is used for 
#' constructing a mip/lp problem.
#' @param E normalized \code{editmatrix}
#' @param prefix \code{character} used for naming dummy variables in matrix.
#' @keywords internal
softEdits.editmatrix <- function(E, prefix="delta.", postfix="", M=1e7, ...){
  
  if (!nrow(E)){
    return(E)
  }
  
  n <- nrow(E)
  vars <- getVars(E)
  ops <- getOps(E)
  
  adapt <- paste(prefix,rownames(E), sep="")
  
  A <- getA(E)
  b <- getb(E)
  isna <- is.na(b)
  eq <- (ops == "==") & !isna
  
  Ab <- cbind( A
             , diag(-M, n)
             , b
             )[!isna,,drop=FALSE]
  
  # copy the equality constraints
  Ab_eq <- if(any(eq)){
           cbind( -A
                , diag(-M,n)
                , -b
                )[eq,,drop=FALSE]
           }
  
  # clear A, trick that keeps the rownames
  A[,] <- 0
  
  # NA's must be changed.
  Ab_na <- if(any(isna)){
              cbind( A #matrix(0, nrow=n, ncol=ncol(A))
                , diag(1, n)
                , 1
                )[isna,,drop=FALSE]
           }
  # TODO cleanup this code
  #print(list(Ab=Ab, Ab_eq=Ab_eq, Ab_na=Ab_na))
  Ab <- rbind(Ab, Ab_eq, Ab_na)
  rownames(Ab) <- make.unique(paste0(rownames(Ab), postfix), "_" )
  ops <- c(ops[!isna], ops[eq])
  ops <- gsub("==", "<=", ops)
  ops <- c(ops, rep("==", sum(isna)))
  
  colnames(Ab) <- c(getVars(E), adapt, "CONSTANT")
  
  seE <- neweditmatrix(Ab, ops=ops)
  seE
}

#' Derive editmatrix with soft constraints. This is a utility function that is used for 
#' constructing a mip/lp problem.
#' @param E normalized \code{editmatrix}
#' @param prefix \code{character} used for naming dummy variables in matrix.
#' @keywords internal
softEdits.cateditmatrix <- function(E, prefix="delta.", postfix="", ...){
  if (!nrow(E)){
    return(E)
  }
  eq <- getOps(E) == "=="
  b <- getb(E)
  
  dummies <- paste(prefix, rownames(E), sep="")
  
  seA <- diag(ifelse(eq, -1, 1), ncol=length(eq), nrow=length(eq))
  colnames(seA) <- dummies
  seA <- cbind(getA(E), seA)
  rownames(seA) <- paste0(rownames(seA), postfix)
  
  binvars <- sapply(colnames(seA), is.character)
  seE <- as.editmatrix(seA, b, getOps(E), binvars=binvars)
  seE
}

#' Derive editmatrix with soft constraints based on boundaries of variables. This is a utility function that is used for 
#' constructing a mip/lp problem.
#' @param E normalized \code{editmatrix}
#' @param prefix \code{character} used for naming dummy variables in matrix.
#' @keywords internal
softEdits.editarray<- function(E, prefix="delta.", postfix="", ...){
  if (!nrow(E)){
    return(E)
  }
  softEdits.cateditmatrix(cateditmatrix(E), prefix=prefix, postfix=postfix, ...)
}

#quick tests

# E <- editmatrix(expression( x - y < 2
#                           , x + y < 5
#                           , x + y == 3
#                           , z == 1
#                           )
#                )
# 
# # set the z == NA
# E[4,ncol(E)] <- NA
# (se <- softEdits(E))
#softEdits(cateditmatrix(c("if (married) adult", "married %in% c(TRUE,FALSE)","married==TRUE")))
