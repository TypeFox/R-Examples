#' Curve PCA
#'
#' Calculate principal directions of a set of curves
#'
#' @param betamean array (n,T) of mean curve
#' @param mu array (n,T) of mean srvf
#' @param K array (2*T,2*T) covariance matrix
#' @param mode Open ("O") or Closed ("C") curves
#' @param no number of components
#' @param N number of samples on each side of mean
#' @return pd list describing principal directions
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_srvf_align(beta[,,1,1:2],maxit=2) # note: use more shapes, small for speed
#' K = curve_karcher_cov(out$betamean, beta[,,1,1:2])
#' pd = curve_principal_directions(out$betamean, out$q_mu, K)
curve_principal_directions <- function(betamean, mu, K, mode="O", no=3, N=5){
    n = nrow(betamean)
    T1 = ncol(betamean)
    out = svd(K)
    U = out$u
    s = out$d
    V = out$v

    qarray = array(list(), c(no,2*N+1))
    qarray1 = array(list(), N)
    qarray2 = array(list(), N)
    pd = array(list(), c(no, 2*N+1))
    pd1 = array(list(), N)
    pd2 = array(list(), N)

    for (m in 1:no){
        princDir = rbind(U[1:T1,m], U[(T1+1):(2*T1),m])
        v = sqrt(s[m])*princDir
        q1 = mu
        epsilon = 2./N

        # forward direction from mean
        for (i in 1:N){
            normv = sqrt(innerprod_q2(v,v))

            if (normv < 1e-4){
                q2 = mu
            } else {
                q2 = cos(epsilon*normv)*q1 + sin(epsilon*normv)*v/normv
                if (mode=="C"){
                    q2 = project_curve(q2)
                }
            }

            qarray1[[i]] = q2
            p = q_to_curve(q2)
            centroid1 = -1 * calculatecentroid(p)
            dim(centroid1) = c(length(centroid1),1)
            out = scale_curve(p+repmat(centroid1,1,T1))
            pd1[[i]] = out$beta_scaled

            # Parallel translate tangent vector
            basis2 = find_basis_normal(q2)
            v = parallel_translate(v, q1, q2, basis2, mode)

            q1 = q2
        }

        # Backward direction from mean
        v = -sqrt(s[m])*princDir
        q1 = mu
        for (i in 1:N){
            normv = sqrt(innerprod_q2(v,v))

            if (normv < 1e-4){
                q2 = mu
            } else {
                q2 = cos(epsilon*normv)*q1 + sin(epsilon*normv)*v/normv
                if (mode=="C"){
                    q2 = project_curve(q2)
                }
            }

            qarray2[[i]] = q2
            p = q_to_curve(q2)
            centroid1 = -1*calculatecentroid(p)
            dim(centroid1) = c(length(centroid1),1)
            out = scale_curve(p+repmat(centroid1, 1, T1))
            pd2[[i]] = out$beta_scaled

            # parallel translate tangent vector
            basis2 = find_basis_normal(q2)
            v = parallel_translate(v, q1, q2, basis2, mode)

            q1 = q2
        }

        for (i in 1:N){
            qarray[[m,i]] = qarray2[[N+1-i]]
            pd[[m,i]] = pd2[[N+1-i]]
        }

        qarray[[m,N+1]] = mu
        centroid1 = -1 * calculatecentroid(betamean)
        dim(centroid1) = c(length(centroid1),1)
        out = scale_curve(betamean + repmat(centroid1,1,T1))

        pd[[m,N]] = out$beta_scaled

        for (i in 1:N){
            qarray[[m,i+N+1]] = qarray1[[i]]
            pd[[m,i+N+1]] = pd1[[i]]
        }
    }

    return(pd)
}
