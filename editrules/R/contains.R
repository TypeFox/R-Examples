#' Determine which edits contain which variable(s)
#'
#' For an \code{\link{editmatrix}}, variables with coefficients smaller than
#' \code{tol} are considered not to be contained in an edit.
#'
#' @param E \code{\link{editarray}}, \code{\link{editmatrix}}, \code{\link{editset}}, or \code{matrix}
#' @param var \code{character}, names of a categorical variables in \code{E}. If var=NULL, all variables are treated.
#' @param ... arguments to be passed to other methods
#' @return \code{logical} vector of length nrow(E), TRUE for edits containing \code{var}
#' @export
contains <- function(E,var=NULL,...){
  if (nedits(E)==0){
    if ( is.null(var) ) var <- getVars(E)
    return(
      array(
        logical(0),
        dim=c(0,length(var)),
        dimnames=list(edit=NULL,variable=var)
      )
    )
  }
  UseMethod('contains')
}


#' matrix method for contains
#' 
#' @method contains matrix
#' @rdname contains
#' @export
#' @keywords internal
contains.matrix <- function(E, var=NULL, tol=sqrt(.Machine$double.eps),...){
    if ( is.null(var) && !is.null(colnames(E)) ) var <- colnames(E)
    if ( is.null(var) ) var <- 1:ncol(E)
    u <- abs(E[,var,drop=FALSE]) > tol
   
    dimnames(u) <- dimnames(E[,var,drop=FALSE])
    u
}


#' editmatrix method for contains
#' 
#' @method contains editmatrix
#' @rdname contains
#' @param tol tolerance to check zero-entries
#' @export
#' @keywords internal
contains.editmatrix <- function(E, var=NULL, tol=sqrt(.Machine$double.eps), ...){
    A <- getA(E)
    vars <- getVars.editmatrix(E)
    if (is.null(var)){ 
      var <- vars
    } else {
      stopifnot(all(var %in% vars))
    }
    if ( length(var) > 0 ){ 
      u <- abs(A[,var,drop=FALSE]) > tol
    } else {
      u <- array(FALSE,dim=c(nrow(A),0))
    } 
    dimnames(u) <- list(edit=rownames(E),variable=var)
    u
}

#' contains method for editarray
#'
#' @method contains editarray
#' @rdname contains
#' @export
#' @keywords internal
contains.editarray <- function(E,var=NULL,...){
    if ( !is.editarray(E) ) stop("Argument not of class editarray")
    ind <- getInd(E)
    if ( is.null(var)) var <- names(ind)
   
    contains.boolmat(getArr(E),ind,var)
}

#' Determine if a boolean matrix contains \code{var}
#'
#' @param A array
#' @param ind index
#' @param var variable name
#' @keywords internal
contains.boolmat <- function(A, ind, var){
    ind <- ind[var]
    v <- vapply(ind, function(ii) rowSums(A[,ii,drop=FALSE]) < length(ii), FUN.VALUE=logical(nrow(A)))
    if ( is.vector(v) ){ 
        v <- array(v,dim=c(1,length(v)), dimnames=list(edit=rownames(A),var=names(v)))
    } else {
        dimnames(v) <- list(edit=rownames(A),variable=colnames(v))
    } 
  v
}

#' contains method for cateditmatrix
#'
#' @method contains editset
#' @rdname contains 
#' @export
#' @keywords internal
contains.cateditmatrix <- function(E, var=NULL, ...){
    m <- contains.editmatrix(E)
    ind <- indFromArray(m,sep=":")
    vapply(ind, function(ii) rowSums(m[,ii,drop=FALSE])>0, FUN.VALUE=logical(nrow(m)) )
}

#' contains method for editset
#'
#' @method contains editset
#' @rdname contains 
#' @export
#' @keywords internal
contains.editset <- function(E,var=NULL,...){

  if ( is.null(var) ) var <- getVars(E)
  nedits <- nrow(E$num) + nrow(E$mixcat)
  # create a boolean array holding the answer
  T <- array(FALSE,
    dim=c(nedits,length(var)),
      dimnames=list(
        edits=c(
          rownames(E$num),
          rownames(E$mixcat)
        ),
        variables=var
      )
  )
  # contains for numerical variables
  numvar <- var[var %in% getVars(E$num)]
  nnum <- nrow(E$num)
  ivr <- var %in% numvar
  if ( any(ivr) && length(numvar) > 0 && nnum > 0 ){  
    T[1:nrow(E$num),numvar] <- contains(E$num, var[ivr])
  }
  # contains for categorical variables in conditional edits (mix)
  nmix <- nrow(E$mixcat)
  if (nmix>0){
    imix <- (nnum+1):(nnum+nmix) 
    catvar <- var[var %in% getVars(E,type='cat')]
    if (length(catvar) > 0){
      T[imix,catvar] <- contains(E$mixcat, catvar) 
    }
    # contains for numerical variables in mixed edits
    vnm <- var[var %in% getVars(E$mixnum)]
    if ( length(vnm) > 0 ){
      vmc <- getVars(E$mixcat)
      emn <- rownames(E$mixnum)
      vmc <- vmc[vmc %in% emn]
      vmc <- vmc[match(emn,vmc)]
      X <- contains(E$mixcat)[,vmc,drop=FALSE]
      Y <- contains(E$mixnum,vnm)
      T[imix,vnm] <- (X%*%Y) > 0
    }
  }
  T
}


