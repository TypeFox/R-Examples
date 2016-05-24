rppvalue <- function(x, method = c("bootstrap", "normal"), x.boot) {
    stopifnot(is.numeric(x), is.character(method))

    method <- match.arg(method)
    if(method != "bootstrap" && !missing(x.boot)){
        warning(paste("The argument x.boot is ignored for method = ", method, ".", sep = ""))
    }
    switch(method,
        normal    = pnorm(x),
        bootstrap = {
            stopifnot(!missing(x.boot))
            stopifnot(is.numeric(x.boot))
            if(!is.matrix(x)){
                x <- matrix(x, nrow = 1L)
            }
            if(!is.matrix(x.boot)){
                x.boot <- matrix(x.boot, nrow = 1L)
            }
            x <- array(x, dim = c(dim(x), dim(x.boot)[2L]))
            return(drop(rowMeans(sweep(x, 2L:3L, x.boot, ">="), dims = 2L)))
        }
    )
}
