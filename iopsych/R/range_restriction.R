# This file contains range restriction formulas


#' Lawley multivariate range restriction correction.
#'
#' @param rcov The covariance matrix of the restricted sample.
#' @param vnp The covariance matrix of predictors explicitly used for selection. 
#'        This matrix should be based on the the unrestricted population.
#' @param as_cor This argument can be set to FALSE to return a covariance matrix.
#' @return The the correlation matrix or variance covariance in the unrestricted population.
#' @author The original function was written by Adam Beatty and adapted by Allen Goebl.
#' @references Lawley D. N (1943). A note on Karl Pearson's selection formulae. 
#'             \emph{Proceedings of the Royal Society of Edinburgh.},
#'             62(Section A, Pt. 1), 28-30.
#' @examples
#' data(rcea1994)
#' vstar <- rcea1994$vstar
#' vpp   <- rcea1994$vpp
#'
#'lMvrrc(rcov=vstar, vnp=vpp)
#' @export
#' @importFrom stats cov2cor 
lMvrrc <- function (rcov, vnp, as_cor=TRUE) {
    #Define index
    pp <- 1:ncol(vnp)
    np <- (max(pp)+1):ncol(rcov)
    #Define submatrices
    pp_pp <- rcov[pp,pp]
    pp_np <- rcov[pp,np]
    np_pp <- rcov[np,pp]
    np_np <- rcov[np,np]
    #Invert pp_pp
    i_pp_pp <- qr.solve(pp_pp)
    #Define output
    Vpp_pp <- vnp
    Vpp_np <- vnp %*% i_pp_pp %*% pp_np
    Vnp_pp <- t(Vpp_np)
    Vnp_np <- np_np + (np_pp %*% i_pp_pp %*% (Vpp_np - pp_np))
    #Return output
    out <- (cbind(rbind(Vpp_pp, Vnp_pp), rbind(Vpp_np, Vnp_np)))
    if(as_cor) { return(cov2cor(out)) }
    return(out)
}

