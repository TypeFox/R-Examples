
checkGH <- function(object) {

    ### check gradient  and hessian
    suppressWarnings(gr <- numDeriv::grad(object$loglik, coef(object), 
                                          weights = weights(object)))
    s <- Gradient(object)
    cat("Compare gradients")
    print(all.equal(gr, s, check.attributes = FALSE))

    suppressWarnings(H1 <- numDeriv::hessian(object$loglik, coef(object), 
                                             weights = weights(object)))
    H2 <- Hessian(object)
    cat("Compare hessians:")
    print(all.equal(H1, H2, check.attributes = FALSE))
}
