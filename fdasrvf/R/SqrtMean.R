#' SRVF transform of warping functions
#'
#' This function calculates the srvf of warping functions with corresponding
#' shooting vectors
#'
#' @param gam matrix (\eqn{N} x \eqn{M}) of \eqn{M} warping functions with \eqn{N} samples
#' @return Returns a list containing \item{mu}{Karcher mean psi function}
#' \item{gam_mu}{Karcher mean warping function}
#' \item{psi}{srvf of warping functions}
#' \item{vec}{shooting vectors}
#' @keywords srvf alignment
#' @references Srivastava, A., Wu, W., Kurtek, S., Klassen, E., Marron, J. S.,
#'  May 2011. Registration of functional data using fisher-rao metric,
#'  arXiv:1103.3817v2 [math.ST].
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_warp")
#' out = SqrtMean(simu_warp$gam)
SqrtMean <- function(gam){
    TT = nrow(gam)
    n = ncol(gam)

    psi = matrix(0,TT-1,n)
    for (i in 1:n){
        psi[,i] = sqrt(diff(gam[,i])*TT)
    }

    # Find Direction
    mnpsi = rowMeans(psi)
    dqq = sqrt(colSums((psi - matrix(mnpsi,ncol=n,nrow=TT-1))^2))
    min_ind = which.min(dqq)
    mu = psi[,min_ind]
    tt = 1
    maxiter = 10
    vec = matrix(0,TT-1,n)
    lvm = rep(0,5)
    for (iter in 1:maxiter){
        for (i in 1:n){
            v = psi[,i] - mu
            dot<- simpson(seq(0,1,length.out=TT-1),mu*psi[,i])
            dot.limited<- ifelse(dot>1, 1, ifelse(dot<(-1), -1, dot))
            len = acos(dot.limited)
            if (len > 0.0001){
                vec[,i] = (len/sin(len))*(psi[,i] - cos(len)*mu)
            }else{
                vec[,i] = rep(0,TT-1)
            }
        }
        vm = rowMeans(vec)
        lvm[iter] = sqrt(sum(vm*vm)/TT)
        mu = cos(tt*lvm[iter])*mu + (sin(tt*lvm[iter])/lvm[iter])*vm
        if (lvm[iter] < 1e-6 || iter >=maxiter){
            break
        }
    }

    phi = matrix(0,TT,n)
    for (i in 1:n){
        tmp = rep(0,TT)
        tmp[1:TT-1] = psi[,i]*psi[,i]/TT
        phi[,i] = cumsum(tmp)
    }
    gam_mu = gam_mu = c(0,cumsum(mu*mu))/n
    gam_mu = (gam_mu - min(gam_mu))/(max(gam_mu)-min(gam_mu))

    out = list(mu = mu,gam_mu = gam_mu,psi = psi,vec = vec)
    return(out)

}
