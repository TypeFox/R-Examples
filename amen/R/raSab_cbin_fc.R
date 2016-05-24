#' Simulate a and Sab from full conditional distributions under the cbin
#' likelihood
#' 
#' Simulate a and Sab from full conditional distributions under the cbin
#' likelihood
#' 
#' 
#' @usage raSab_cbin_fc(Z, Y, a, b, Sab, odmax, odobs, SS =
#' round(sqrt(nrow(Z))))
#' @param Z a square matrix, the current value of Z
#' @param Y square matrix of ranked nomination data
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
#' @export raSab_cbin_fc
raSab_cbin_fc<-
function (Z, Y, a, b, Sab, odmax, odobs, SS = round(sqrt(nrow(Z))))
{
    E <- Z - a %*% t(rep(1, nrow(Z)))
    MEL <- MEU <- -E
    MEL[!is.na(Y) & Y == 0] <- -Inf
    MEU[!is.na(Y) & Y == 1] <- Inf
    MEL[is.na(Y)] <- -Inf
    MEU[is.na(Y)] <- Inf
    lba <- apply(MEL, 1, max)
    lba[is.na(lba)] <- -Inf
    uba <- apply(MEU, 1, min)
    uba[is.na(uba)] <- Inf
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

