translifeExp <- function(coef, var, ns, pfixed){
    ## Just 'beta <-- -beta'
    namn <- names(coef)
    if (pfixed){
        A <- diag(c(rep(-1, length(coef) - ns), rep(1, ns)))
    }else{
        A <- diag(c(rep(-1, length(coef) - 2 * ns), rep(1, 2 * ns)))
    }
    coef <- as.vector(A %*% coef)
    var <- A %*% var %*% t(A)

    names(coef) <- namn
    colnames(var) <- rownames(var) <- namn
    
    list(coefficients = coef, var = var)
}
