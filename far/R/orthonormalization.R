#  *****************************************************************************
#   File : orthonormalization.R
#         ************************************************************
#   Description : 
#       Orthonormalization of a set of vector
#       using the Gram-Schmidt algorithm
#   Version : 1.1
#   Date : 2003-04-06
#         ************************************************************
#   Author : Julien Damon <julien.damon@gmail.com>
#   License : LGPL
#   URL: https://github.com/Looping027/far
#  *****************************************************************************

# ******************************************************************************
#   Title : orthonormalization
#         ************************************************************
#   Description : Orthonormalization of a set of vector
#       using the Gram-Schmidt algorithm
#   Version : 1.1
#   Date : 2003-04-06
# ******************************************************************************
orthonormalization <- function(u=NULL,basis=TRUE,norm=TRUE)
{
    if (is.null(u)) return(NULL)
    # change the type of u if necessary
    if (!(is.matrix(u)) )
        u <- as.matrix(u)
    p <- nrow(u)  # dimension of the space
    n <- ncol(u)  # number of vectors
    if (prod(abs(La.svd(u)$d)>1e-8)==0)
        stop("colinears vectors")
    if (p < n)  # test the number of vectors
    {
        warning("too much vectors to orthogonalize.")
        u <- as.matrix(u[,1:p])
        n <- p
    }

    if (basis & (p > n)) # if a basis is desired
                         # u is extended to be square
    {
        base <- diag(p) # the canonical basis of size p x p
        coef.proj <- crossprod(u,base) / diag(crossprod(u))
        # calculation of base2,
        # the orthogonal (to the span space of u) part of the base
        base2 <- base - u %*% matrix(coef.proj,nrow=n,ncol=p)
        # norms of the vector of base2
        norm.base2 <- diag(crossprod(base2))
        base <- as.matrix(base[,order(norm.base2) > n])
        # extention of u with the more orthogonal vector of base
        u <- cbind(u,base)
        n <- p
    }

    v <- u  # initialization
    # start of the gram-schmidt algorithm
    if (n > 1)
    {
        for (i in 2:n)
        {
            coef.proj <- c(crossprod(u[,i],v[,1:(i-1)])) /
                    diag(crossprod(v[,1:(i-1)]))
            v[,i] <- u[,i] - matrix(v[,1:(i-1)],nrow=p) %*% 
                                matrix(coef.proj,nrow=i-1)
        }
    }
    if (norm) # if orthonormal vector are wanted
    {
        coef.proj <- 1/sqrt(diag(crossprod(v)))
        v <- t(t(v) * coef.proj)
    }
    return(v)
}
