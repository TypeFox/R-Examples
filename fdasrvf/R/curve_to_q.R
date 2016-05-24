#' Convert to SRVF space
#'
#' This function converts curves to SRVF
#'
#' @param beta array describing curve (n,T)
#' @return q array describing srvf
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' q = curve_to_q(beta[,,1,1])
curve_to_q <- function(beta){
    n = nrow(beta)
    T1 = ncol(beta)
    v = matrix(0, n, T1)
    for (i in 1:n){
        v[i,] = gradient(beta[i,], 1.0/(T1-1))
    }

#     len = sum(sqrt(colSums(v*v)))/T1
#     v = v/len
    q = matrix(0,n,T1)
    for (i in 1:T1){
        L = sqrt(pvecnorm(v[,i],2))
        if (L>0.0001){
            q[,i] = v[,i]/L
        } else {
            q[,i] = v[,i]*0.0001
        }
    }

    q = q/sqrt(innerprod_q2(q, q))

    return(q)
}
