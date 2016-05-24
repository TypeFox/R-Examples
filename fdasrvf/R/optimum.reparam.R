#' Align two functions
#'
#' This function aligns two SRSF functions using Dynamic Programming
#'
#' @param Q1 srsf of function 1
#' @param T1 sample points of function 1
#' @param Q2 srsf of function 2
#' @param T2 sample points of function 2
#' @param lambda controls amount of warping (default = 0)
#' @param method controls which optmization method (default="DP") options are
#' Dynamic Programming ("DP"), Coordinate Descent ("DP2"), Riemannian BFGS
#' ("RBFGS") and Simultaneous Alignment ("SIMUL")
#' @param w controls LRBFGS (default = 0.01)
#' @param f1o initial value of f1, vector or scalar depending on q1, defaults to zero
#' @param f2o initial value of f2, vector or scalar depending on q1, defaults to zero
#' @return gam warping function
#' @keywords srsf alignment, pca
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' q = f_to_srvf(simu_data$f,simu_data$time)
#' gam = optimum.reparam(q[,1],simu_data$time,q[,2],simu_data$time)
optimum.reparam <- function(Q1,T1,Q2,T2,lambda=0,method="DP",w=0.01,f1o=0.0,
                            f2o=0.0){
    n = length(T1)
    Q1=(Q1/pvecnorm(Q1,2))
    Q2=(Q2/pvecnorm(Q2,2))
    C1=srsf_to_f(Q1,T1,f1o)
    C2=srsf_to_f(Q2,T2,f2o)
    rotated = FALSE
    isclosed = FALSE
    skipm = 0
    auto = 0
    if (method=="DP"){
        G = rep(0,n)
        T = rep(0,n)
        size = 0;
        ret = .Call('DPQ2', PACKAGE = 'fdasrvf', Q1, T1, Q2, T2, 1, n, n, T1, T2, n, n, G, T, size, lambda);

        G = ret$G[1:ret$size]
        Tf = ret$T[1:ret$size]
        gam0 = approx(Tf,G,xout=T2)$y
    } else if (method=="SIMUL"){
        out = simul_align(C1,C2)
        u = seq(0,1,length.out=length(out$g1))
        tmin = min(T1)
        tmax = max(T1)
        timet2 = T1
        timet2 = (timet2-tmin)/(tmax-tmin)
        gam0 = simul_gam(u,out$g1,out$g2,timet2,out$s1,out$s2,timet2)
    } else if (method=="DP2") {
        opt = rep(0,n+1+1);
        swap = FALSE
        fopts = rep(0,5)
        comtime = rep(0,5)

        out = .Call('opt_reparam', PACKAGE = 'fdasrvf', C1,C2,n,1,0.0,TRUE,
                    rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)
        gam0 = out$opt
        gam0 = gam0[1:(length(gam0)-2)]

        if (out$swap){
            gam0 = invertGamma(gam0);
        }
    } else {
        opt = rep(0,n+1+1);
        swap = FALSE
        fopts = rep(0,5)
        comtime = rep(0,5)

        out = .Call('opt_reparam', PACKAGE = 'fdasrvf', C1,C2,n,1,w,FALSE,
                    rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)

        if (out$fopts[1] == 1000){
            out = .Call('opt_reparam', PACKAGE = 'fdasrvf', C1,C2,n,1,0.0,TRUE,
                        rotated,isclosed,skipm,auto,opt,swap,fopts,comtime)
        }

        gam0 = out$opt
        gam0 = gam0[1:(length(gam0)-2)]

        if (out$swap){
            gam0 = invertGamma(gam0);
        }
    }

    gam = (gam0-gam0[1])/(gam0[length(gam0)]-gam0[1])  # slight change on scale

    return(gam)
}
