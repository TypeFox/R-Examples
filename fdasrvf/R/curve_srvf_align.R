#' Align Curves
#'
#' Aligns a collection of curves using the elastic square-root velocity (srvf) framework.
#'
#' @param beta array (n,T,N) for N number of curves
#' @param mode Open ("O") or Closed ("C") curves
#' @param maxit maximum number of iterations
#' @return Returns a list containing \item{betan}{aligned curves}
#' \item{qn}{aligned srvfs}
#' \item{betamean}{mean curve}
#' \item{q_mu}{mean SRVFs}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_srvf_align(beta[,,1,1:2],maxit=2) # note: use more shapes, small for speed
#' K = curve_karcher_cov(out$betamean, beta[,,1,1:2])
curve_srvf_align <- function(beta, mode="O", maxit=20){
    tmp = dim(beta)
    n = tmp[1]
    T1 = tmp[2]
    N = tmp[3]
    out = curve_karcher_mean(beta, mode, maxit)
    mu = out$mu
    betamean = out$betamean
    v = out$v
    q = out$q

    qn = array(0, c(n,T1,N))
    betan = array(0, c(n,T1,N))
    centroid2 = calculatecentroid(betamean)
    dim(centroid2) = c(length(centroid2),1)
    betamean = betamean - repmat(centroid2, 1, T1)
    q_mu = curve_to_q(betamean)

    # align to mean
    for (ii in 1:N){
        beta1 = beta[,,ii]
        centroid1 = calculatecentroid(beta1)
        dim(centroid1) = c(length(centroid1),1)
        beta1 = beta1 - repmat(centroid1,1,T1)

        # Iteratively optimize over SO(n) x Gamma
        # Optimize over SO(n) x Gamma
        out = reparam_curve(betamean, beta1)
        gamI = invertGamma(out$gam)
        beta1 = out$R %*% shift_f(beta1, out$tau)

        # Apply optimal re-parameterization to the second curve
        beta1 = group_action_by_gamma_coord(beta1, gamI)

        # Optimize over SO(n)
        out = find_rotation_seed_coord(betamean, beta1)
        qn[,,ii] = curve_to_q(out$beta2new)
        betan[,,ii] = out$beta2new
    }

    return(list(betan=betan, qn=qn, betamean=betamean, q_mu=q_mu))
}
