## Auxiliary coefficients
alpha <- function(model) {
    KL1     <- K(model)
    adjcoef <- uniroot(f         = KL1,
                       interval  = c(.Machine$double.eps^0.75,
                                     min(model$rates) - .Machine$double.eps),
                       difforder = 0L,
                       tol       = .Machine$double.eps)$root

    bound <- 1.0 / KL1(adjcoef, 1L)

    ## y = horizon / x
    function(y) {
        if (y > bound) {
            return(c(par   = -adjcoef,
                     value = 0.0))
        } else {
            alpha <- optim(par    = 0.0,
                           fn     = function(arg) KL1(arg, 0L) - arg / y,
                           gr     = function(arg) KL1(arg, 1L) - 1.0 / y,
                           method = 'BFGS')[['par']]
            return(c(par   = -alpha,
                     value = KL1(alpha)))
        }
    }
}
