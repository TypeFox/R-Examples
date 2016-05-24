contr.XuWu <- function (n, contrasts=TRUE) 
{
    ## the contrasts option does not do anything but is needed for model.matrix
    ## to work on objects of this type
    if (!contrasts) stop("contr.XuWu not defined for contrasts=FALSE")

    ## function to calculate orthogonal normalized contrasts 
    ## that satisfy the XuWu normalization
    ##      Xu and Wu call these orthonormal
    ##      however, this is confusing
    ## These are based on Helmert contrasts
    if (length(n) <= 1L) {
        if (is.numeric(n) && length(n) == 1L && n > 1L) 
            levels <- seq_len(n)
        else stop("not enough degrees of freedom to define contrasts")
    }
    else levels <- n
    levels <- as.character(levels)
        n <- length(levels)
        cont <- array(-1, c(n, n - 1L), list(levels, NULL))
        for (j in 1:(n-1))  
        cont[1:j, j] <- - sqrt(n/(j*(j+1)))
        cont[col(cont) <= row(cont) - 2L] <- 0
        cont[col(cont) == row(cont) - 1L] <- sqrt(n*seq_len(n - 1L)/(1+seq_len(n - 1L)))
        colnames(cont) <- NULL
   cont         
}

contr.XuWuPoly <- function (n, contrasts=TRUE) 
{
    ## the contrasts option does not do anything but is needed for model.matrix
    ## to work on objects of this type
    
    if (!contrasts) stop("contr.XuWuPoly not defined for contrasts=FALSE")

    ## function to calculate orthogonal normalized contrasts 
    ## that satisfy the XuWu normalization
    ##      Xu and Wu call these orthonormal
    ##      however, this is confusing
    ## here based on polynomial contrasts
    cont <- contr.poly(n)
    cont * sqrt(nrow(cont))
}
