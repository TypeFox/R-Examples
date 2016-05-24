saddlepointVaru <- function(process, type = 2L) {

    switch(type,
        { # TYPE 1 (switching back and forth between the monetary domain and the frequency
          # domain)
            vx  <- process[['vx']]
            phi <- function(x, prob, i) {
                if (i == 0L) {
                    return(qexp(prob, adjcoef(process), lower.tail = FALSE))
                } else {
                    v <- vx(x)
                    return(x + 0.5 * (qnorm(prob)^2.0 - v$z^2.0) / v$v)
                }
            }
            varu <- function(prob, n = 4L) {
                stopifnot(n >= 0L)
                for(i in 0L:n) {
                    x <- phi(x, prob, i)
                }
                return(x)
            }
        },
        { # TYPE 2 (iteration in the frequency domain)
            zv    <- process[['zv']]
            KL.d1 <- process[['KL.d1']]
            KL.d2 <- process[['KL.d2']]

            phi <- function(x, prob, i) {
                if (i == 0L) {
                    return(adjcoef(process) * (1.0 + 1.0 / log(prob)))
                } else {
                    if (i == 1L) {
                        div <- qexp(prob, adjcoef(process), lower.tail = FALSE)
                    } else {
                        div <- x * KL.d2(x)
                    }
                    return(x + 0.5 * (qnorm(prob)^2.0 - zv(x)^2.0) / div)
                }
            }
            varu <- function(prob, n = 4L) {
                stopifnot(n >= 0L)
                for(i in 0L:n) {
                    x <- phi(x, prob, i)
                }
                return(structure(.Data       = KL.d1(x),
                                 saddlepoint = x))
            }
        }
    )

    return(varu)
}
