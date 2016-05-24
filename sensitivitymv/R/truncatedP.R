truncatedP <-
function (p, trunc = 0.2) 
    {
        stopifnot((trunc>0)&(trunc<=1))
        stopifnot(is.vector(p))
        stopifnot((min(p)>=0)&(max(p)<=1))
        w <- prod(p^(p<=trunc))
        L <- length(p)
        if (w > trunc) {
            1
        }
        else {
            pr <- 0
            for (k in 1:L) {
                s <- 0:(k - 1)
                term1 <- sum(w * (w <= (trunc^k)) * (((k * log(trunc)) - 
                                                          log(w))^s)/factorial(s))
                term2 <- (trunc^k) * (w > (trunc^k))
                pr <- pr + (choose(L, k) * ((1 - trunc)^(L - k)) * 
                                (term1 + term2))
            }
            pr
        }
    }
