saddlepointTvaru <- function(process, type = c('tail', 'density'), ...) {
    type <- match.arg(type)

    psi    <- attr(saddlepointRuinprob(process, ...), 'diagnostics')
    myvaru <- saddlepointVaru(process, type = 2L)
    myadj  <- adjcoef(process)

    function(prob, n = 4L) {
        quant  <- myvaru(prob, n)
        intmin <- attr(quant, 'saddlepoint')
        quant  <- as.vector(quant)
        p.len  <- length(prob)
        aux    <- lower.tri(x    = matrix(data = NA_real_,
                                          ncol = p.len,
                                          nrow = p.len),
                            diag = TRUE)

### FIX ME: use crossprod() or tcrossprod() instead of %*%
        switch(type,
               tail = {
### FIX ME: psi$psi.v() is not returned anymore by saddlepointRuinprob()
                   drop(int.multi(f = function(v) {
                                      psi$psi.v(v) * process[['KL.d2']](v)
                                  },
                                  nodes = c(intmin, myadj)) %*% aux) / prob + quant
               },
               density = {
### FIX ME: psi$psi.1.v() is not returned anymore by saddlepointRuinprob()
                   drop(int.multi(f = function(v) {
                                      psi$psi.1.v(v) * process[['KL.d1']](v) * process[['KL.d2']](v)
                                  },
                                  nodes = c(intmin, myadj)) %*% aux * process[['p']] * process[['zeta']]) / prob
               })
    }
}
