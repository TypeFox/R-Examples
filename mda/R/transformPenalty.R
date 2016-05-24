transformPenalty <-
function (Q, prior, cl, df = NULL, tot.df = NULL) 
{
    if (missing(Q)) 
        Q <- meanPenalty(prior, cl)
    if (missing(prior)) 
        prior <- attr(Q, "prior")
    if (missing(cl)) 
        cl <- attr(Q, "cl")
    transform.pen <- function(Q, prior, df) {
        df.inv <- function(d, df, lambda = NULL, iterations = 10) {
            if (is.null(lambda)) {
                lambda <- 0.1
                while (sum(1/(1 + d * lambda)) >= df) lambda <- lambda * 
                  2
            }
            df.diriv <- function(d, lambda) -sum((d * lambda)/(1 + 
                d * lambda)^2)
            current.df <- sum(1/(1 + d * lambda))
            if (abs((df - current.df)/df) < 1e-04 | iterations == 
                1) 
                return(list(lambda = lambda, df = current.df))
            else {
                lambda <- exp(log(lambda) - (current.df - df)/df.diriv(d, 
                  lambda))
                Recall(d, df, lambda, iterations - 1)
            }
        }
        pQp <- Q/outer(sqrt(prior), sqrt(prior))
        d <- svd(pQp)$d
        lambda <- df.inv(d, df)$lambda
        lambda * Q
    }
    if (!is.null(tot.df)) {
        if (tot.df >= length(prior)) 
            return(Q * 0)
        else return(transform.pen(Q, prior, tot.df))
    }
    else {
        ncl <- unique(cl)
        df <- rep(df, length = length(ncl))
        for (i in seq(along = ncl)) {
            which <- cl == ncl[i]
            Q[which, which] <- Recall(Q[which, which, drop = FALSE], 
                prior[which], tot.df = df[i])
        }
        return(Q)
    }
}

