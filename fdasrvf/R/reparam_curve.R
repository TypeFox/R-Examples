#' Align two curves
#'
#' This function aligns two SRVF functions using Dynamic Programming
#'
#' @param beta1 array defining curve 1
#' @param beta2 array defining curve 1
#' @param lambda controls amount of warping (default = 0)
#' @param method controls which optmization method (default="DP") options are
#' Dynamic Programming ("DP"), Coordinate Descent ("DP2"), Riemannian BFGS
#' ("RBFGS")
#' @param w controls LRBFGS (default = 0.01)
#' @param rotated boolean if rotation is desireid
#' @param isclosed boolean if curve is closed
#' @return return a List containing \item{gam}{warping function}
#' \item{R}{rotation matrix}
#' \item{tau}{seed point}
#' @keywords srvf alignment
#' @references Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.
#' @export
#' @examples
#' data("mpeg7")
#' gam = reparam_curve(beta[,,1,1],beta[,,1,5])
reparam_curve <- function(beta1,beta2,lambda=0,method="DP",w=0.01,rotated=T,
                          isclosed=F){

    n1 = nrow(beta2)
    M = ncol(beta2)
    timet = seq(0,1,length.out=M)
    skipm = 4
    auto = 2
    tau = 0
    if (method=="DP"){
        # Optimze over SO(n) x Gamma
        q1 = curve_to_q(beta1)

        # Optimzie over SO(n)
        out = find_rotation_seed_coord(beta1, beta2);
        beta2 = out$beta2
        R = out$O_hat
        tau = out$tau
        q2 = curve_to_q(beta2)

        # Optimzie over Gamma
        q1i = q1
        dim(q1i) = c(M*n1)
        q2i = q2
        dim(q2i) = c(M*n1)
        G = rep(0,M)
        T1 = rep(0,M)
        size = 0;
        ret = .Call('DPQ2', PACKAGE = 'fdasrvf', q1i, timet, q2i, timet, n1, M, M, timet, timet, M, M, G, T1, size, lambda);

        G = ret$G[1:ret$size]
        Tf = ret$T[1:ret$size]
        gam0 = approx(Tf,G,xout=timet)$y
    } else if (method=="DP2") {
        c1 = t(beta1)
        dim(c1) = c(M*n1)
        c2 = t(beta2)
        dim(c2) = c(M*n1)
        opt = rep(0,M+n1*n1+1);
        swap = FALSE
        fopts = rep(0,5)
        comtime = rep(0,5)

        out = .Call('opt_reparam', PACKAGE = 'fdasrvf', c1,c2,M,n1,0.0,TRUE,
                    rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)

        tmp = length(out$opt)
        gam0 = out$opt[1:(tmp-5)]
        R = matrix(out$opt[(tmp-4):(tmp-1)],nrow=2)

        if (out$swap){
            gam0 = invertGamma(gam0);
            R = t(R)
        }

    } else {
        c1 = t(beta1)
        dim(c1) = c(M*n1)
        c2 = t(beta2)
        dim(c2) = c(M*n1)
        opt = rep(0,M+n1*n1+1);
        swap = FALSE
        fopts = rep(0,5)
        comtime = rep(0,5)

        out = .Call('opt_reparam', PACKAGE = 'fdasrvf', c1,c2,M,n1,w,FALSE,
                    rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)

        if (out$fopts[1] == 1000){
            out = .Call('opt_reparam', PACKAGE = 'fdasrvf', c1,c2,M,n1,0.0,TRUE,
                        rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)
        }

        tmp = length(out$opt)
        gam0 = out$opt[1:(tmp-5)]
        R = matrix(out$opt[(tmp-4):(tmp-1)],nrow=2)

        if (out$swap){
            gam0 = invertGamma(gam0);
            R = t(R)
        }
    }

    gam = (gam0-gam0[1])/(gam0[length(gam0)]-gam0[1])  # slight change on scale

    return(list(gam=gam,R=R,tau=tau))
}
