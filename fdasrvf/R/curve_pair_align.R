#' Pairwise align two curves
#'
#' This function aligns to curves using Elastic Framework
#'
#' @param beta1 array describing curve 1 (n,T)
#' @param beta2 array describing curve 2 (n,T)
#' @return a list containing \item{beta2n}{aligned curve 2 to 1}
#' \item{q2n}{aligned srvf 2 to 1}
#' \item{gam}{warping function}
#' \item{q1}{srvf of curve 1}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_pair_align(beta[,,1,1], beta[,,1,5])
curve_pair_align <- function(beta1, beta2){
    T1 = ncol(beta1)
    centroid1 = calculatecentroid(beta1)
    dim(centroid1) = c(length(centroid1),1)
    beta1 = beta1 - repmat(centroid1, 1, T1)
    centroid2 = calculatecentroid(beta2)
    dim(centroid2) = c(length(centroid2),1)
    beta2 = beta2 - repmat(centroid2, 1, T1)

    q1 = curve_to_q(beta1)

    # Iteratively optimize over SO(n) x Gamma using old DP
    out = reparam_curve(beta1, beta2)
    beta2n = out$R %*% shift_f(beta2, out$tau)
    gamI = invertGamma(out$gam)
    beta2n = group_action_by_gamma_coord(beta2n, gamI)
    out = find_rotation_seed_coord(beta1, beta2n)
    q2n = curve_to_q(out$beta2new)

    return(list(beta2n=out$beta2new, q2n=q2n, gam=gamI, q1=q1))
}
