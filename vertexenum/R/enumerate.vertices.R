##' @useDynLib vertexenum
##' @export
enumerate.vertices <- function(A, b, warn_if_open=FALSE) {
    APPROXIMATION = 10000;
    ## Given a linear system Ax <= b, return the vertices
    ## of the system

    ## If the matrix or vector contains real numbers, we need to approximate
    ## them with rationals!
    ## Now build the input file for lrs
    m <- dim(A)[1]
    n <- dim(A)[2]
    if (m < 2) {
        print("Matrix must have at least 2 row")
        return(NA)
    }
    if (n < 2) {
        print("Matrix must have at least 2 cols")
        return(NA)
    }
    A_num <- matrix(, nrow=m, ncol=n)
    A_den <- matrix(, nrow=m, ncol=n)
    b_num <- vector(mode="integer", length=m)
    b_den <- vector(mode="integer", length=m)
    for (i in 1:m) {
        ## First approximate b
        cur = numbers::ratFarey(b[i], APPROXIMATION, upper=FALSE)
        b_num[i] = cur[1]
        b_den[i] = cur[2]
        ## Now approximate A
        for (j in 1:n) {
            cur = numbers::ratFarey(A[i,j], APPROXIMATION, upper=FALSE)
            A_num[i,j] = -cur[1]
            A_den[i,j] = cur[2]
        }
    }
    ## print(A_num)
    ## print(A_den)
    ## print(b_num)
    ## print(b_den)
    ## call the external function
    output <- .Call("vertexenum", as.integer(A_num), as.integer(A_den), as.integer(b_num), as.integer(b_den), dim(A))
    ## prune out extreme rays, if they exist
    if (!is.na(match(0, output[,1])) && warn_if_open) {
        print("Input system contains an extreme ray.")
    }
    output <- output[output[,1] != 0,]
    if (is.null(dim(output))) {
        return(output[2:(n+1)])
    } else {
        return(output[,2:(n+1)])
    }
}
