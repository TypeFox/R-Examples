#' Sample shapes from model
#'
#' @param mu array (n,T) of mean srvf
#' @param K array (2*T,2*T) covariance matrix
#' @param mode Open ("O") or Closed ("C") curves
#' @param no number of principal components
#' @param numSamp number of samples
#' @return samples list of sample curves
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_srvf_align(beta[,,1,1:2],maxit=2) # note: use more shapes, small for speed
#' K = curve_karcher_cov(out$betamean, beta[,,1,1:2])
#' samples = sample_shapes(out$q_mu, K)
sample_shapes <- function(mu, K, mode="O", no=3, numSamp=10){
    n = nrow(mu)
    T1 = ncol(mu)

    out = svd(K)
    U = out$u
    s = out$d
    V = out$v

    if (mode == "O"){
        N = 2
    } else {
        N = 10
    }

    epsilon = 1./(N-1)

    q1 = mu
    q2 = mu
    samples = array(list(),numSamp)

    for (i in 1:numSamp){
        v = matrix(0, 2, T1)
        for (m in 1:no){
            v = v + rnorm(1)*sqrt(s[m])*c(U[1:T1,m], U[(T1+1):(2*T1),m])
        }

        q1 = mu
        for (j in 1:(N-1)){
            normv = sqrt(innerprod_q2(v,v))

            if (normv < 1e-4){
                q2 = mu
            } else {
                q2 = cos(epsilon*normv)*q1+sin(epsilon*normv)*v/normv
                if (mode == "C"){
                    q2 = project_curve(q2)
                }
            }

            # Parallel translate tanent vector
            basis2 = find_basis_normal(q2)
            v = parallel_translate(v, q1, q2, basis2, mode)

            q1 = q2
        }

        samples[[i]] = q_to_curve(q2)
    }

    return(samples)
}
