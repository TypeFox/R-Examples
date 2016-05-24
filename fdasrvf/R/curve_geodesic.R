#' Form geodesic between two curvese
#'
#' Form geodesic between two curves using Elastic Method
#'
#' @param beta1 array describing curve 1 (n,T)
#' @param beta2 array describing curve 2 (n,T)
#' @param k number of curves along geodesic (default 5)
#' @return a list containing \item{geod}{curves along geodesic (n,T,k)}
#' \item{geod_q}{srvf's along geodesic}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_geodesic(beta[,,1,1], beta[,,1,5])
curve_geodesic <- function(beta1, beta2, k=5){
    n = nrow(beta1)
    T1 = ncol(beta1)
    beta1 = resamplecurve(beta1, T1)
    beta2 = resamplecurve(beta2, T1)
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

    # form geodesic between the registered curves
    dist = acos(innerprod_q2(q1, q2n))
    geod = array(0, c(n,T1,k))
    geod_q = array(0, c(n,T1,k))

    for (tau in 1:k){
        s = dist*(tau-1)/(k-1)
        geod_q[,,tau] = (sin(dist-s)*q1+sin(s)*q2n)/sin(dist)
        geod[,,tau] = q_to_curve(geod_q[,,tau])
    }

    return(list(geod=geod,geod_q=geod_q))
}
