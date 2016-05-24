readFullMat <- function(file, ...) {
    my.df <- read.table(file = file, ...)
    no.lines <- nrow(my.df)
    no.var <- ncol(my.df)
    if (no.lines%%no.var == 0) {
        no.groups <- no.lines/no.var
    } else {
        stop("No. of lines read is not divided by the no. of variables.")
    }
    var.names <- paste("x", 1:no.var, sep = "")
    my.list <- split(my.df, rep(1:no.groups, each = no.var))
#    my.mat <- lapply(my.list, matrix, byrow = TRUE, ncol = no.var, nrow = no.var)
    # Add variable names into the matrices
    out <- lapply(my.list, function(x, v.names) {
        dimnames(x) <- list(v.names, v.names)
        x
    }, var.names)
    out
    lapply(out, function(x) {as.matrix(x)} )
}



readLowTriMat <- function(file, no.var, ...) {
    if (missing(no.var)) 
        stop("No. of variables was missing.")
    # problem: read by row major!
    my.scan <- scan(file = file, ...)
    ps <- no.var * (no.var + 1)/2
    no.groups <- length(my.scan)/ps
    if (length(my.scan)%%ps != 0) 
        stop("No. of elements read != no.var*(no.var+1)*no.of.studies.")
    my.df <- matrix(my.scan, ncol = ps, nrow = no.groups, byrow = TRUE)
    my.list <- split(my.df, 1:no.groups)
    # mat 1: no.var by no.var of 0
    # mat 2: fill lower triangle of mat by my.list
    # mat 3: lower + upper
    my.mat <- lapply(my.list, function(x, no.var) {
        mat <- matrix(0, ncol = no.var, nrow = no.var)
        mat[upper.tri(mat, diag = TRUE)] <- x
        mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
        mat
    }, no.var)
    # Add variable names into the matrices
    var.names <- paste("x", 1:no.var, sep = "")
    out <- lapply(my.mat, function(x, v.names) {
        dimnames(x) <- list(v.names, v.names)
        x
    }, var.names)
    out
}


#### Convert list of matrices into row vector of matrices
#my.mat <- do.call(rbind, lapply(my.df, vech))
#write.table(my.mat, file='row.dat', row.names=FALSE, col.names=FALSE)
# http://tolstoy.newcastle.edu.au/R/help/04/11/6694.html

#### Read cov/cor elements by column major
#### Sample data of 2 studies with 3 variables
#### Study1: c11 c21 c31 c22 c32 c33
#### Study2: c11 c21 c31 c22 c32 c33
readStackVec <- function(file, ...) {
    my.df <- as.matrix(read.table(file = file, ...))
    # Convert no. of covariance elements into no. of variables by p*=p(p+1)/2
    no.col <- ncol(my.df)
    no.var <- (sqrt(1 + 8 * no.col) - 1)/2
    if (abs(no.var - round(no.var)) > .Machine$double.eps^0.5) {
        stop("No. of columns does not match no. of variables: p(p+1)/2")
    }
    no.groups <- nrow(my.df)
    my.list <- split(my.df, 1:no.groups)

    # mat 1: no.var by no.var of 0
    # mat 2: fill lower triangle of mat by my.list
    # mat 3: lower + upper
    my.mat <- lapply(my.list, function(x, no.var) {
        mat <- matrix(0, ncol = no.var, nrow = no.var)
        mat[lower.tri(mat, diag = TRUE)] <- x
        mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
        mat
        }, no.var)
    
    # Add variable names into the matrices
    var.names <- paste("x", 1:no.var, sep = "")
    out <- lapply(my.mat, function(x, v.names) {
        dimnames(x) <- list(v.names, v.names)
        x
    }, var.names)
    out
}
