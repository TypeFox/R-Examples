#' Creates a feasible starting matrix with a desired mean average degree
#'
#' This extension of \code{\link{findFeasibleMatrix}} attempts to
#' create a feasible matrix where a certain proportion of the entries
#' is positive. There is no guarantee that this proportion is
#' achieved. If it is not possible then this matrix will report a warning
#' and simply return the matrix constructed by findFeasibleMatrix.
#'
#' @inheritParams findFeasibleMatrix
#' @param targetmean Average proportion of positive entries of the resulting matrix. Defaults to 0.3
#' @export
findFeasibleMatrix_targetmean <- function(r,c,p,eps=1e-9,targetmean=0.3){
    Lorig <- findFeasibleMatrix(r=r,c=c,p=p,eps=eps)
    if (mean(Lorig>0)<targetmean){
        avgavailable <- mean(p>0)
        psample <- min(1,(targetmean-mean(Lorig>0))/avgavailable)
        Adjstart <- matrix(rbinom(length(r)*length(c),
                                   size=1,
                                   prob=psample*(p>0)),
                         nrow=length(r),
                         ncol=length(c)
                         )
        maxsaveval <- pmin(matrix(r/pmax(1,rowSums(Adjstart)),nrow=length(r),ncol=length(c)),
                           matrix(c/pmax(1,colSums(Adjstart)),nrow=length(r),ncol=length(c),byrow=TRUE))
        Lstart <- matrix(runif(length(r)*length(c))*maxsaveval,nrow=length(r))*Adjstart
        Lstart+findFeasibleMatrix(r=r-rowSums(Lstart),c=c-colSums(Lstart),p=p,eps=eps)
    }else{
        warning("Desired mean degree is less than minimal degree that is necessary.");
        Lorig
    }
}


#' Creates a feasible starting matrix
#'
#' Creates a matrix with nonnegative entries, given row and column
#' sums and 0 on the diagonal. Superseeded by the more flexible findFeasibleMatrix.
#'
#' @param L Vector of row sums
#' @param A Vector of column sums
#' @return A matrix with nonnegative entries and given row/column sums
#' and 0 on the diagonal.
#' @examples
#' getfeasibleMatr(c(0.5,1,0),c(0.5,0,1))
#' getfeasibleMatr(rep(1,4),rep(1,4))
#' getfeasibleMatr(2^(1:3),2^(3:1))
#' getfeasibleMatr(1:5,1:5)
#' getfeasibleMatr(1:5,5:1)
#' @export
getfeasibleMatr <- function(L,A){
    if (min(L,A)<0) stop("L and A must be nonnegative.")
    if (abs(sum(L)-sum(A))>1e-10) stop("sum(L) different to  sum(A).")
    if (max(A+L)>sum(L)) stop("No feasible matrix exists: max(A+L)>sum(L)")
    Lorig <- L; Aorig <- A
    n <- length(L)
    resall <- matrix(0,nrow=n,ncol=n)
    while(sum(L)>1e-10&&sum(A)>1e-10){
        if (min(ifelse(A>0,A,Inf))<min(ifelse(L>0,L,Inf))){
            transp <- TRUE
            Atmp <- A
            A <- L
            L <- Atmp
        }else{
            transp <- FALSE
        }
        i <- which.min(ifelse(L>0,L,Inf))
        LA <- ifelse(A>0,L+A,0)
        LA[i] <- Inf
        o <- order(LA,decreasing=TRUE)
        L <- L[o]
        A <- A[o]
        LA <- LA[o]
        xi=max(which((LA)[-1]==(LA)[2]))
        if (xi+1==n)
            Delta <- min(L[1],xi*(LA)[2])
        else
            Delta <- min(L[1],xi*((LA)[2]-(LA)[2+xi]))
        res <- matrix(0,nrow=n,ncol=n)
        res[1,2:(1+xi)]=Delta/xi
        L[1] <- L[1]-Delta
        A[2:(1+xi)] <- A[2:(1+xi)]-Delta/xi
        res[o,o] <- res
        if (transp)res <- t(res)
        resall <- resall+res
        A[o] <- A ##sort back
        L[o] <- L ##sort back
        if (transp){
            Atmp <- A
            A <- L
            L <- Atmp
        }
    }
    if (any(abs(rowSums(resall)-Lorig)>1e-9)||any(abs(colSums(resall)-Aorig)>1e-9)){
        warning("Discrepancy between intended row/column sums and ",
                "row/column sums of generated matrix detected. ",
                "This may only be a numerical issue. ",
                "Maximum absolute difference between row and column sums: ",
                max(abs(rowSums(resall)-Lorig),abs(colSums(resall)-Aorig)))
    }
    resall
}
