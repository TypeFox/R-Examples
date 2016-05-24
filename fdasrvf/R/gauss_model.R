#' Gaussian model of functional data
#'
#' This function models the functional data using a Gaussian model extracted from
#' the principal components of the srvfs
#'
#' @param fn matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param qn matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned srvfs
#' @param gam warping functions
#' @param n number of random samples (n = 1)
#' @param sort_samples sort samples (default = F)
#' @return Returns a list containing \item{fs}{random aligned samples}
#' \item{gams}{random warping function samples}
#' \item{ft}{random function samples}
#' @keywords pca
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_data")
#' data("simu_warp")
#' out1 = gauss_model(simu_warp$fn,simu_data$time,simu_warp$qn,simu_warp$gam,n = 10)
gauss_model <- function(fn,time,qn,gam,n = 1,sort_samples = FALSE){
    # Parameters
    no = 3
    eps = .Machine$double.eps
    binsize = mean(diff(time))
    M = length(time)

    # compute mean and covariance in q-domain
    mq_new = rowMeans(qn)
    m_new = sign(fn[round(length(time)/2),])*sqrt(abs(fn[round(length(time)/2),]))  # scaled version
    mqn = c(mq_new,mean(m_new))
    C = cov(t(rbind(qn,m_new)))

    q_s = rmvnorm(n,mean=mqn,sigma=C,method="svd")
    q_s = t(q_s)
    end = dim(q_s)[1]

    # compute the correspondence to the original function domain
    fs = matrix(0,M,n)
    for (k in 1:n){
        fs[,k] = cumtrapzmid(time,q_s[1:(end-1),k]*abs(q_s[1:(end-1),k]),sign(q_s[end,k])*(q_s[end,k]^2))
    }

    # random warping generation
    rgam = randomGamma(gam,n)
    gams = matrix(0,n,M)
    for (k in 1:n){
        gams[k,] = invertGamma(rgam[k,])
    }
    gams = t(gams)

    # sort functions and warpings
    if (sort_samples == T){
        mx = apply(fs,2, max)
        out_sort = sort(mx,index.return=TRUE)
        seq1 = out_sort$ix

        # compute the psi-function
        fy = gradient(t(rgam),binsize)
        psi = fy/sqrt(abs(fy)+eps)
        psi = t(psi)
        ip = rep(0,n)
        len = rep(0,n)
        for (i in 1:n){
            ip[i] = rep(1,M)%*%psi[i,]/M;
            len[i] = acos(rep(1,M)%*%psi[i,]/M)
        }
        out_sort = sort(len,index.return=TRUE)
        seq2 = out_sort$ix

        # combine x-variability and y-variability
        ft = matrix(0,M,n)
        for (k in 1:n){
            tmp = approx((0:(M-1))/(M-1),fs[,seq1[k]],xout = gams[,seq2[k]])
            ft[,k] = tmp$y
            while (is.na(ft[,k])){
                rgam2 = randomGamma(gam,1)
                tmp = approx((0:(M-1))/(M-1),fs[,seq1[k]],xout = invertGamma(rgam2))
                ft[,k] = tmp$y
            }
        }
    }else
    {
        # combine x-variability and y-variability
        ft = matrix(0,M,n)
        for (k in 1:n){
            tmp = approx((0:(M-1))/(M-1),fs[,k],xout = gams[,k])
            ft[,k] = tmp$y
            while (is.na(ft[,k])[1]){
                rgam2 = randomGamma(gam,1)
                tmp = approx((0:(M-1))/(M-1),fs[,k],xout = invertGamma(rgam2))
                ft[,k] = tmp$y
            }
        }
    }

    samples = list(fs = fs, gams = rgam, ft = ft)

    return(samples)
}
