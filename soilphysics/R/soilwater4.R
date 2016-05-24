soilwater4 <- 
function(psi, Bd, a, b, c, model = c("Silva", "Ross"))
{
    model <- match.arg(model)
    if (model == "Silva") {
       if (length(psi) != length(Bd)) 
           stop("incompatible dimensions!")
       dat <- cbind(psi, Bd)
       if (!is.numeric(dat)) 
           stop("non-numeric data!")
       theta <- exp(a + b * Bd) * psi^c
    } else {
       stopifnot(inherits(psi, "numeric"))
       theta <- a * psi^c
    }
    return(theta)
}