##
##  p ol y a d d .R  Polynomial Addition
##


polyadd <- function(p, q){
    if ( (!is.vector(p, mode="numeric") && !is.vector(p, mode="complex")) ||
         (!is.vector(q, mode="numeric") && !is.vector(q, mode="complex")) )
        stop("Arguments 'p' and 'q' must be real or complex vectors.")

    lp <- length(p)
    lq <- length(q)

    if (lp >= lq) {
        r <- p + c(numeric(lp-lq), q)
    } else {
        r <- q + c(numeric(lq-lp), p)
    }

    lr <- length(r)
    while (r[1] == 0 && lr > 1) {
        r <- r[2:lr]
        lr <- lr - 1
    }
    return(r)
}
