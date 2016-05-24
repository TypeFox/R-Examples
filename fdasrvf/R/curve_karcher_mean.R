#' Karcher Mean of Curves
#'
#' Calculates Karcher mean of a collection of curves using the elastic square-root velocity (srvf) framework.
#'
#' @param beta array (n,T,N) for N number of curves
#' @param mode Open ("O") or Closed ("C") curves
#' @param maxit maximum number of iterations
#' @return Returns a list containing \item{mu}{mean srvf}
#' \item{betamean}{mean curve}
#' \item{v}{shooting vectors}
#' \item{q}{array of srvfs}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' out = curve_karcher_mean(beta[,,1,1:2],maxit=2) # note: use more shapes, small for speed
curve_karcher_mean <- function(beta, mode="O", maxit=20){
    tmp = dim(beta)
    n = tmp[1]
    T1 = tmp[2]
    N = tmp[3]
    q = array(0, c(n,T1,N))
    for (ii in 1:N){
        q[,,ii] = curve_to_q(beta[,,ii])
    }

    # Initialize mu as one of the shapes
    mnq = rowMeans(q[1,,])
    dqq = sqrt(colSums((q[1,,] - matrix(mnq,ncol=N,nrow=T1))^2))
    min_ind = which.min(dqq)
    mu = q[,,min_ind]
    betamean = beta[,,min_ind]

    delta = 0.5
    tolv = 1e-4
    told = 5*1e-3
    itr = 1
    sumd = rep(0,maxit+1)
    v = array(0,c(n,T1,N))
    normvbar = rep(0,maxit+1)

    while (itr<maxit){
        cat(sprintf("Iteration: %d\n",itr))

        mu = mu / sqrt(innerprod_q2(mu,mu))

        sumv = matrix(0,2,T1)
        sumd[itr] = 0.

        # TODO: parallelize
        for (i in 1:N){
            out = karcher_calc(beta[,,i], q[,,i], betamean, mu, mode)
            v[,,i] = out$v
            sumd[itr+1] = sumd[itr+1] + out$d^2
        }

        sumv = rowSums(v,dims=2)

        # compute average direction of tangent vectors v_i
        vbar = sumv/N

        normvbar[itr] = sqrt(innerprod_q2(vbar,vbar))
        normv = normvbar[itr]

        if ((normv>tolv) && (abs(sumd[itr+1]-sumd[itr])>told)){
            # update mu in direction of vbar
            mu = cos(delta*normvbar[itr])*mu + sin(delta*normvbar[itr])*vbar/normvbar[itr]

            if (mode=="C"){
                mu = project_curve(mu)
            }

            x = q_to_curve(mu)
            a = -1 * calculatecentroid(x)
            dim(a) = c(length(a),1)
            betamean = x + repmat(a,1,T1)
        } else {
            break
        }

        itr = itr + 1
    }

    return(list(mu=mu,betamean=betamean,v=v,q=q))
}
