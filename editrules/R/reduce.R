
#' Remove redundant variables and edits.
#'
#' Remove variables which are not contained in any edit and remove edits which are 
#' \code{\link[=isObviouslyRedundant]{obviously redundant}}.
#'
#' @param E \code{\link{editmatrix}} or \code{\link{editarray}}
#' @param ... arguments to pass to other methods
#'
#' @export
#' @seealso \code{\link{contains}}, \code{\link{eliminate}}, \code{\link{substValue}}
#'
reduce <- function(E,...){
    UseMethod('reduce')
}


#' 
#' @method reduce editmatrix
#' @param tol elements of \code{E} with absolute value < \code{tol} are considered 0.
#' 
#' @rdname reduce
#' @export
reduce.editmatrix  <- function(E, tol=sqrt(.Machine$double.eps),...){
  if ( nrow(E) == 0 ) return(neweditmatrix(matrix(numeric(0)),character(0)))
  m <- as.matrix(E)
  if ( tol > 0 ) m[abs(m) < tol] <- 0
  B <- m != 0
  v <- 1:(ncol(m)-1)
  vars <- which(colSums(B[,v,drop=FALSE]) != 0)
  edits <- (rowSums(B) != 0)
  E[edits,c(vars,ncol(m)) , drop=FALSE]
}


#'
#' @method reduce editarray
#' 
#' @export
#' @rdname reduce
reduce.editarray <- function(E,...){
    E <- E[!isObviouslyRedundant.editarray(E),,drop=FALSE]
    m <- as.matrix(E)
    ind <- getInd(E)
    H <- getH(E)
    b <- sapply(ind,function(ind) all(m[,ind,drop=FALSE])) 
    if ( any(b) ){
        J <- logical(0)   
        for ( j in ind[b] ) J <- c(J,j)
        sep=getSep(E)
        m <- m[,-J,drop=FALSE]
        ind <- indFromArray(m,sep=sep)
        i <- apply(!m,1,all)
        m <- m[!i,,drop=FALSE]
        if (!is.null(H))  H <- H[!i,,drop=FALSE]
        E <- neweditarray(E=m, ind=ind, sep=sep, names=rownames(m), H=H)
    }
    E
}


#'
#' @method reduce editset
#' @export
#' @rdname reduce
#'
reduce.editset <- function(E,...){
    num <- reduce(E$num)
    mixcat <- reduce(E$mixcat)
    v <- getVars(mixcat)
    mixnum <- reduce(E$mixnum[rownames(E$mixnum) %in% v,])

    imix <- grepl("^.l",v)
    if ( nrow(mixcat) > 0 ){
        m <- logical(nedits(mixcat)) 
        if ( any(imix) ) m <- apply(contains(mixcat, v[imix]), 1, any)
        pref <- ifelse(m,"mix","cat")
        rownames(mixcat) <- paste(pref, 1:nrow(mixcat), sep="")
    }
    simplify(neweditset(
        num = num,
        mixnum = mixnum,
        mixcat = mixcat
    ))

}

#' Retrieve values stricktly implied by rules
#'
#' @export
impliedValues <- function(E,...){
  UseMethod('impliedValues')
}

#'
#'
#' Detects cases where two inequalities imply an equality, e.g. \eqn{x\leq 0} and \eqn{x\geq0}
#' implies \eqn{x=0}. Also detects straight equalities, e.g. \eqn{x==0} implies \eqn{x=0}. Such
#' cases arise frequently when manipulating edits by value subsitution or variable elimination.
#' The function recursively detects equalities and combined inequalities that imply fixed values, 
#' substitutes those fixed values and looks for new implied values until no new values are found.
#'
#' @method impliedValues editmatrix
#'
#' @param E editmatrix
#' @param tol Maximum deviation for two values to be considered equal.
#' @param ... Currently unused
#' @return Numeric vector, whose names are variable names and values are unique values implied by the rules.
#' @rdname impliedValues
#'
#' @seealso \code{\link{reduce}}, \code{\link{substValue}}, \code{\link{eliminate}}
#' 
#' @export
impliedValues.editmatrix <- function(E,tol=sqrt(.Machine$double.eps),...){
  if (!isNormalized(E)) E <- normalize(E)
 
  iv <- numeric(0)
  iv[getVars(E)] <- NA
 
  # edits containing one variable
  i <- rowSums(contains(E,tol=tol)) == 1
  if ( sum(i) == 0 ) return(numeric(0)) 
  e <- reduce(E[i,])
 
  # direct equalities
  ops <- getOps(e)
  eqs <- e[ops == '==', ]
  eqA <- getA(eqs)
  eqb <- getb(eqs)
  j <- apply(abs(eqA),1,which.max)
  i <- 1:nrow(eqs)
  iv[colnames(eqA)[j]] <- eqb/eqA[cbind(i,j)]
 
  x <- iv[!is.na(iv)]
  # combined inequalities
  ieq <- e[ops %in% c('<=','>=')]
  if (length(x) > 0)
    ieq <- substValue(ieq, names(x), x,reduce=TRUE)
  Ab <- getAb(e)
  for ( v in getVars(ieq) ){
    ab <- Ab[abs(Ab[,v])>tol,,drop=FALSE]
    ab <- ab/abs(ab[,v])
    b <- ab[,'CONSTANT']
    for ( i in 1:length(b) ){
      j <- ( abs(b[i] + b[-i]) < tol & abs(ab[i,v] + ab[-i,v]) < tol )
      if (any(j)){
        iv[v] = abs(b[i])
        break
      }
    }
  }
  iv <- iv[!is.na(iv)]
  if ( length(iv) > 0 ){
    E <- substValue(E,names(iv),iv,reduce=TRUE)
    c(iv,impliedValues.editmatrix(E,tol))
  } else {
    iv
  }
}
