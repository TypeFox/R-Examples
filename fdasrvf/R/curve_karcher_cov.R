#' Curve Karcher Covariance
#'
#' Calculate Karcher Covariance of a set of curves
#'
#' @param betamean array (n,T) of mean curve
#' @param beta array (n,T,N) for N number of curves
#' @param mode Open ("O") or Closed ("C") curves
#' @return K covariance matrix
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_srvf_align(beta[,,1,1:2],maxit=2) # note: use more shapes, small for speed
#' K = curve_karcher_cov(out$betamean, beta[,,1,1:2])
curve_karcher_cov <- function(betamean, beta, mode="O"){
    tmp = dim(beta)
    n = tmp[1]
    T1 = tmp[2]
    N = tmp[3]

    # Compute Karcher covariance of uniformly sampled mean
    betamean = resamplecurve(betamean, T1)
    mu = curve_to_q(betamean)
    if (mode=="C"){
        mu = project_curve(mu)
        basis = find_basis_normal(mu)
    }

    v = array(0, c(n,T1,N))
    for (i in 1:N){
        beta1 = beta[,,i]

        out = inverse_exp_coord(betamean, beta1)
        # Project to the tangent space of manifold to obtain v_i
        if (mode=="O"){
            v[,,i] = out$v
        } else {
            v[,,i] = project_tangent(out$v, mu, basis)
        }
    }

    K = matrix(0, 2*T1, 2*T1)

    for (i in 1:N){
        w = v[,,i]
        w = c(w[1,], w[2,])
        K = K + w %*% t(w)
    }

    K = K / (N-1)

    return(K)
}
