phypoexp <- function(q, rate = 1.0, lower.tail = TRUE, log.p = FALSE, tailarea = FALSE) {

    if (tailarea) {
        mycoef <- ratetoalpha(rate) / (rate * sum(1.0 / rate))
    } else {
        mycoef <- ratetoalpha(rate)
    }

    res <- drop(tcrossprod(outer(X          = q,
                                 Y          = rate,
                                 FUN        = pexp,
                                 lower.tail = lower.tail,
                                 log.p      = FALSE),
                           t(mycoef)))

    if (log.p) {
        return(log(res))
    } else {
        return(res)
    }
}
