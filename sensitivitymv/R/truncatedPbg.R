truncatedPbg <-
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
            mix <- dbinom(1:L,L,trunc)
            prob <- 1 - pgamma(-log(w/(trunc^(1:L))),1:L)
            sum(mix*prob)
        }
    }
