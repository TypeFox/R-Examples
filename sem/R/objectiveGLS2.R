# last modified 2012-01-06 by J. Fox

objectiveGLS2 <- function (gradient = FALSE) 
{
    result <- list(objective = function(par, model.description) {
        with(model.description, {
            A <- P <- matrix(0, m, m)
            val <- ifelse(fixed, ram[, 5], par[sel.free])
            A[arrows.1] <- val[one.head]
            P[arrows.2t] <- P[arrows.2] <- val[!one.head]
            I.Ainv <- solve(diag(m) - A)
            C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
#            Cinv <- solve(C)
            SS <- invS %*% (S - C)
            f <- 0.5 * sum(diag(SS %*% SS))
            attributes(f) <- list(C = C, A = A, P = P)
            f
        })
    })
    class(result) <- "semObjective"
    result
}