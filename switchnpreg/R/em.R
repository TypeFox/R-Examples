switchnpreg <- function(x, y, f, alpha, sigma2, lambda, ...,
                        method = c('pl', 'bayes'),
                        var.equal = TRUE,
                        z.indep = TRUE,
                        eps.cv, eps.em, maxit.cv, maxit.em) {
    method <- match.arg(method)
    
    fj_update <- switch(method,
                        bayes = fj_bayes,
                        pl = fj_pl)
    
    sigma_update <- if (var.equal) equal.sigma2.update else diff.sigma2.update

    alpha_update <- if (z.indep) alpha_iid else alpha_markov
    pij_update <-  if (z.indep) pij_iid else pij_markov
    criteria_alpha <- if (z.indep) criteria_alpha_iid else criteria_alpha_markov
    standard_error <- if (z.indep) stderr_iid else stderr_markov
    
    em(x, y, f, alpha, sigma2, lambda, ...,
       eps.cv = eps.cv,
       eps.em = eps.em,
       maxit.cv = maxit.cv,
       maxit.em = maxit.em,
       pijFn = pij_update,
       updateFjFn = fj_update,
       updateAlphaFn = alpha_update,
       updateSigmaFn = sigma_update,

       criteriaAlphaFn = criteria_alpha,
       stderrFn = standard_error)
}


em <- function(x, y, f, alpha, sigma2,
               lambda, ...,
               eps.cv, eps.em,
               maxit.cv, maxit.em,
               pijFn,
               updateFjFn,
               updateAlphaFn,
               updateSigmaFn,
               criteriaAlphaFn,
               stderrFn) {
    n <- length(y)

    J <- length(alpha)
    
    c.cv  <- rep(Inf, J)

    current <- list(f = f,
                    sigma2 = sigma2,
                    alpha = alpha,
                    pij = NULL)

    for (iter.cv in seq(length = maxit.cv)) {
        c.em <- rep(Inf, length(eps.em))

        for (iter.em in seq(length = maxit.em)) {
            updates <- emStep(x = x, y = y,
                              current = current,
                              lambda = lambda,
                              ...,
                              pijFn = pijFn,
                              updateFjFn = updateFjFn,
                              updateAlphaFn = updateAlphaFn,
                              updateSigmaFn = updateSigmaFn)

            c.em <- compare(current, updates,
                            criteria.alpha.function = criteriaAlphaFn)

            ## updating
            current <- updates

            if (all(c.em <= eps.em)) break
        }

        c.cv <- colMeans(abs(f - current$f))

        if (any(c.cv > eps.cv)) {
            lambda <- updateLambda(updateFjFn,
                                   x = x,
                                   y = y,
                                   current = current,
                                   ...)
        } else break
        
        f <- current$f
    }

    stderr <- stderrFn(y, current)
    
    list(current = current,
         lambda = lambda,
         iter.cv = iter.cv,
         stderr = stderr)
}


### Performs a single EM step
emStep <- function(x, y, current, lambda,
                   ...,
                   pijFn,
                   updateFjFn,
                   updateAlphaFn,
                   updateSigmaFn) {
    
    pij <- pijFn(y, current)
    
    alpha.new <- updateAlphaFn(y, pij, current)
    
    f.new <- updateF(x, y, pij, current, lambda, updateFjFn, ...)
    
    sigma2.new <- updateSigmaFn(y, pij, f.new$f_hat, f.new$traces)
    
    list(f = f.new$f_hat,
         sigma2 = sigma2.new,
         alpha = alpha.new,
         pij = pij)
}
