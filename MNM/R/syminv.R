syminv <- function(X)
        {
        ch.X <- chol(X)
        chol2inv(ch.X)
        }
