"posterior.ante"<-function(x, k=1){

    if (is.null(dim(x))[1]) {
        x <- as.matrix(x)
    }
    if (sqrt(dim(x)[2])%%1 != 0) {
        stop("cannot coerce rows of x into a matrix - the number of columns should be a triangular number")
    }

    if(!k<sqrt(ncol(x))){stop("k must be less than the dimensions of the matrix")}

    coerce.ante <- function(x, k=k) {
        n<-sqrt(length(x))
        V<-matrix(x, n, n)
        cholV<-chol(V)
        sd<-diag(diag(cholV))
        beta<-t(solve(cholV,sd))
        c(diag(sd)^2, unlist(sapply(1:k, function(k){-beta[(1:(n-k)-1)*n+(1:(n-k))+k]})))
    }
    ante<-t(apply(x, 1, coerce.ante, k=k))

    if (is.mcmc(x) == FALSE) {
        warning("posterior.ante expecting mcmc object")
        ante
    }
    else {
        as.mcmc(ante)
    }
}



