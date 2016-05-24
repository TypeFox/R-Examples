hypoexpRuinprob <- function(process) {
    stopifnot(is.hypoexp(process[['claims']]))

    mypoly.factors <- as.polylist(lapply(X   = process[[c('claims', 'hypoexp', 'rates')]],
                                         FUN = function(arg) {
                                             c(arg, -1.0)
                                         }))


    mypoly.rhs <- mean(process[['claims']]) * polynom(c(process[['zeta']], -1.0)) * prod(mypoly.factors)

    mypoly.lhs <- process[['zeta']] * process[['q']] * sum(as.polylist(
        lapply(X   = seq_along(mypoly.factors),
               FUN = function(index) {
                   process[[c('claims', 'hypoexp', 'coef')]][index] * prod(mypoly.factors[-index])
               })
    ))

    r <- solve(mypoly.lhs - mypoly.rhs)

    const <- solve(a = rbind(outer(X   = process[[c('claims', 'hypoexp', 'rates')]],
                                   Y   = r,
                                   FUN = function(.rates, .r) {
                                       .rates / (.rates - .r)
                                   }),
                             rep.int(1.0, length(r))),
                   b = rep.int(1.0, length(r)))

    const1 <- r * const / (process[['p']] * process[['zeta']])
    const2 <- const - const1

    genexp <- function(multarg, exparg, cutoff) {
        function(x) {
            pmin.int(cutoff, Re(drop(crossprod(exp(outer(-exparg, x)), multarg))))
        }
    }

    return(structure(.Data       = list(psi   = genexp(const,  r, 1.0),
                                        psi.1 = genexp(const1, r, 1.0),
                                        psi.2 = genexp(const2, r, 1.0),
                                        dens  = genexp(const * r, r, Inf)),
                     compmethod  = 'hypoexp',
                     riskproc    = process,
                     parameters  = list(NULL),
                     diagnostics = list(C  = const,
                                        C1 = const1,
                                        C2 = const2,
                                        r  = r)))
}
