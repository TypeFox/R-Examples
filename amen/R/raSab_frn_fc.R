#' Simulate a and Sab from full conditional distributions under frn likelihood
#' 
#' Simulate a and Sab from full conditional distributions under frn likelihood
#' 
#' 
#' @usage raSab_frn_fc(Z, Y, YL, a, b, Sab, odmax, odobs, SS =
#' round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square matrix of ranked nomination data
#' @param YL list of ranked individuals, from least to most preferred in each
#' row
#' @param a current value of row effects
#' @param b current value of column effects
#' @param Sab current value of Cov(a,b)
#' @param odmax a scalar or vector giving the maximum number of nominations for
#' each individual
#' @param odobs observed outdegree
#' @param SS number of iterations
#' @return \item{Z}{new value of Z} \item{Sab}{new value of Sab} \item{a}{new
#' value of a}
#' @author Peter Hoff
#' @export raSab_frn_fc
raSab_frn_fc<-
function (Z, Y, YL, a, b, Sab,  odmax, odobs, SS = round(sqrt(nrow(Z))))
{
    E <- Z - a %*% t(rep(1, nrow(Z)))
    lba <- -E[cbind(1:nrow(Z), YL[, 1])]
    lba[is.na(lba)] <- -Inf
    uba <- -apply(E - (Y != 0) * (Inf^(Y != 0)), 1, max, na.rm = TRUE)
    uba[odobs == odmax] <- Inf
    for (ss in 1:SS) {
        ea <- b * Sab[1, 2]/Sab[2, 2]
        sa <- sqrt(Sab[1, 1] - Sab[1, 2]^2/Sab[2, 2])
        a <- ea + sa * qnorm(runif(nrow(Z), pnorm((lba - ea)/sa),
            pnorm((uba - ea)/sa)))
        Sab <- solve(rwish(solve(diag(2) + crossprod(cbind(a,
            b))), 3 + nrow(Z)))
    }
    list(Z = E + a %*% t(rep(1, nrow(Z))), a = a, Sab = Sab)
}

