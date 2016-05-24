#' Bring an (edit) matrix to reduced row echelon form.
#'
#' If \code{E} is a matrix, a matrix in reduced row echelon form is returned.
#' If \code{E} is an \code{\link{editmatrix}} the equality part of \code{E} is transformed
#' to reduced row echelon form. For an \code{\link{editset}}, the numerical part is
#' transformed to reduced row echelon form.
#'
#' @aliases echelon.editmatrix echelon.matrix
#'
#' @param E a matrix or editmatrix
#' @param ... options to pass on to further methods.
#' @export
#' @seealso \code{\link{eliminate}}, \code{\link{substValue}}
echelon <- function(E,...){
    UseMethod("echelon")
}



#' @method echelon editmatrix
#' @rdname echelon
#' @export
echelon.editmatrix <- function(E,...){
    o <- getOps(E)
    # nothing to eliminate?
    eq <- o == "=="
    if ( sum(eq) <= 1 ) return(E)
    Ab <- getAb(E)
    Ab <- rbind(
        echelon.matrix(Ab[eq,,drop=FALSE]),
        Ab[!eq,,drop=FALSE]
    )
    neweditmatrix(Ab,c(o[eq],o[!eq]))

}

#' @rdname echelon
#' @method echelon matrix
#' @param tol tolerance that will be used to determine if a coefficient equals zero.
#' @export
echelon.matrix <- function(E, tol=sqrt(.Machine$double.eps), ...){
    k <- min(ncol(E),nrow(E))
    I <- 1:nrow(E)
    for ( i in 1:k ){
        I1 <- which(I >= i)
        ip <- I1[which.max(abs(E[I1,i]))]
        p <- E[ip,]
        if ( abs(p[i]) < tol ) next
        if ( ip > i ) E[c(ip,i),] <- E[c(i,ip),]
        E[-i,] <- E[-i,] - outer(E[-i,i],p/p[i])
    }
    d <- diag(E)
    id <- abs(d) > tol
    E[id,] <- E[id,]/d[id]
    I0 <- rowSums(abs(E) < tol) == ncol(E)
    rbind(E[!I0,,drop=FALSE],E[I0,,drop=FALSE])
}

#' @method echelon editset
#' @rdname echelon
#' @export
#'
echelon.editset <- function(E,...){
    E$num <- echelon(E$num)
    E
}



