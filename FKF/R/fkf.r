fkf <- function(a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt, check.input = TRUE)
{
    if(any(c(missing(a0), missing(P0), missing(dt), missing(ct), missing(Tt),
             missing(Zt), missing(HHt), missing(GGt), missing(yt)))){

      stop("None of the input arguments 'a0', 'P0', 'dt', 'ct', 'Tt', 'Zt',",
           "'HHt', 'GGt', and 'yt' must be missing.")
    }

    ## 'check.input' should always be 'TRUE' unless the performance
    ## becomes crucial and correctness of the arguments concerning
    ## dimensions, class and storage.mode is ensured.
    if(check.input){

        ## Check the storage mode: Must be 'double' for all arguments
        stor.mode <- c(storage.mode(a0), storage.mode(P0), storage.mode(dt),
                       storage.mode(ct), storage.mode(Tt), storage.mode(Zt),
                       storage.mode(HHt), storage.mode(GGt), storage.mode(yt))

        names(stor.mode) <- c("a0", "P0", "dt", "ct", "Tt", "Zt",
                              "HHt", "GGt", "yt")

        if(any(stor.mode != "double")){
            stop("storage mode of variable(s) '",
                 paste(names(stor.mode)[stor.mode != "double"],
                       collapse = "', '"),
                 "' is not 'double'!\n", sep = "")
        }

        ## Check classes of arguments
        error.string <- ""
        if(!is.vector(a0)){
            error.string <- paste(error.string,
                                  "'a0' must be a vector!\n", sep = "")
        }

        if(!is.matrix(P0)){
            error.string <- paste(error.string,
                                  "'P0' must be of class 'matrix'!\n", sep = "")
        }

        if(!is.matrix(dt) && !is.vector(dt)){
            error.string <- paste(error.string,
                                  "'dt' must be of class 'vector' or 'matrix'!\n", sep = "")
        }else if(is.vector(dt)){
            dt <- as.matrix(dt)
        }

        if(!is.matrix(ct) && !is.vector(ct)){
            error.string <- paste(error.string,
                                  "'ct' must be of class 'vector' or 'matrix'!\n", sep = "")
        }else if(is.vector(ct)){
            ct <- as.matrix(ct)
        }

        if(!is.array(Tt)){
            error.string <- paste(error.string,
                                  "'Tt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(Tt)){
            Tt <- array(Tt, c(dim(Tt), 1))
        }

        if(!is.array(Zt)){
            error.string <- paste(error.string,
                                  "'Zt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(Zt)){
            Zt <- array(Zt, c(dim(Zt), 1))
        }

        if(!is.array(HHt)){
            error.string <- paste(error.string,
                                  "'HHt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(HHt)){
            HHt <- array(HHt, c(dim(HHt), 1))
        }

        if(!is.array(GGt)){
            error.string <- paste(error.string,
                                  "'GGt' must be of class 'matrix' or 'array'!\n", sep = "")
        }else if(is.matrix(GGt)){
            GGt <- array(GGt, c(dim(GGt), 1))
        }

        if(!is.matrix(yt)){
            error.string <- paste(error.string,
                                  "'yt' must be of class 'matrix'!\n", sep = "")
        }
        if(error.string != ""){
            stop(error.string)
        }
        ## Check compatibility of dimensions
        n <- ncol(yt)
        d <- nrow(yt)
        m <- length(a0)

        if(dim(P0)[2] != m | dim(P0)[1] != m | dim(dt)[1] != m |
           dim(Zt)[2] != m | dim(HHt)[1] != m | dim(HHt)[2] != m  |
           dim(Tt)[1] != m  | dim(Tt)[2] != m)
        {
            stop("Some of dim(P0)[2], dim(P0)[1], dim(dt)[1],\n",
                 "dim(Zt)[2], dim(HHt)[1], dim(HHt)[2],\n",
                 "dim(Tt)[1] or dim(Tt)[2] is/are not equal to 'm'!\n")
        }

        if((dim(dt)[2] != n && dim(dt)[2] != 1) |
           (dim(ct)[2] != n && dim(ct)[2] != 1) |
           (dim(Tt)[3] != n && dim(Tt)[3] != 1) |
           (dim(Zt)[3] != n && dim(Zt)[3] != 1) |
           (dim(HHt)[3] != n && dim(HHt)[3] != 1) |
           (dim(GGt)[3] != n && dim(GGt)[3] != 1) |
           dim(yt)[2]  != n)
        {
            stop("Some of dim(dt)[2], dim(ct)[2], dim(Tt)[3],\n",
                 "dim(Zt)[3], dim(HHt)[3], dim(GGt)[3] or\n",
                 "dim(yt)[2] is/are neither equal to 1 nor equal to 'n'!\n")
        }

        if(dim(ct)[1] != d | dim(Zt)[1] != d  |
           dim(GGt)[1]!= d  | dim(GGt)[2] != d  | dim(yt)[1] != d)
        {
            stop("Some of dim(ct)[1], dim(Zt)[1], dim(GGt)[1],\n",
                 "dim(GGt)[2] or dim(yt)[1] is/are not equal to 'd'!\n")

        }

    }

    time.0 <- proc.time()

    ans <-.Call("FKF", a0, P0, dt, ct, Tt, Zt, HHt, GGt, yt, PACKAGE = "FKF")

    ans$sys.time <- proc.time() - time.0

    return(ans)
}

