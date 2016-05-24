#' Vertical Functional Principal Component Analysis
#'
#' This function calculates vertical functional principal component analysis
#' on aligned data
#'
#' @param fn matrix (\eqn{N} x \eqn{M}) of \eqn{M} aligned functions with \eqn{N} samples
#' @param time vector of size \eqn{N} describing the sample points
#' @param qn matrix (\eqn{N} x \eqn{M}) of \eqn{M} of aligned srvfs
#' @param no number of prinicpal components to extract
#' @param showplot show plots of prinipal directions (default = T)
#' @return Returns a list containing \item{q_pca}{srvf principal directions}
#' \item{f_pca}{f principal directions}
#' \item{latent}{latent values}
#' \item{coef}{coefficients}
#' \item{U}{eigenvectors}
#' @keywords srvf alignment
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  May 2012. Generative Models for Function Data using Phase and Amplitude Separation,
#'  submitted to Computational Statistics and Data Analysis.
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_warp")
#' data("simu_data")
#' vfpca = vertFPCA(simu_warp$fn,simu_data$time,simu_warp$qn,no = 3)
vertFPCA <- function(fn,time,qn,no,showplot = TRUE){
    # Parameters
    coef = -2:2
    NP = 1:no  # number of principal components
    Nstd = length(coef)

    # FPCA
    mq_new = rowMeans(qn)
    m_new = sign(fn[round(length(time)/2),])*sqrt(abs(fn[round(length(time)/2),]))  # scaled version
    mqn = c(mq_new,mean(m_new))
    K = cov(t(rbind(qn,m_new))) #out$sigma

    out = svd(K)
    s = out$d
    stdS = sqrt(s)
    U = out$u

    # compute the PCA in the q domain
    q_pca = array(0,dim=c((length(mq_new)+1),Nstd,no))
    for (k in NP){
        for (i in 1:Nstd){
            q_pca[,i,k] = mqn + coef[i]*stdS[k]*U[,k]
        }
    }

    # compute the correspondence to the original function domain
    f_pca = array(0,dim=c((length(mq_new)),Nstd,no))
    for (k in NP){
        for (i in 1:Nstd){
            f_pca[,i,k] = cumtrapzmid(time,q_pca[1:(dim(q_pca)[1]-1),i,k]*
                abs(q_pca[1:(dim(q_pca)[1]-1),i,k]),sign(q_pca[dim(q_pca)[1],i,k])*
                (q_pca[dim(q_pca)[1],i,k]^2))
        }
    }
    N2 = dim(qn)[2]
    c = matrix(0,N2,no)
    for (k in NP){
        for (i in 1:N2){
            c[i,k] = sum((c(qn[,i],m_new[i])-mqn)*U[,k])
        }
    }

    vfpca = list()
    vfpca$q_pca = q_pca
    vfpca$f_pca = f_pca
    vfpca$latent = s
    vfpca$coef = c
    vfpca$U = U

    if (showplot){
        layout(matrix(c(1,2,3,4,5,6), 2, 3, byrow = TRUE))
        dims = dim(q_pca)
        matplot(time,q_pca[1:(dims[1]-1),,1],type="l")
        title(main="q domain: PD 1")
        matplot(time,q_pca[1:(dims[1]-1),,2],type="l")
        title(main="q domain: PD 2")
        matplot(time,q_pca[1:(dims[1]-1),,3],type="l")
        title(main="q domain: PD 3")
        matplot(time,f_pca[,,1],type="l")
        title(main="f domain: PD 1")
        matplot(time,f_pca[,,2],type="l")
        title(main="f domain: PD 2")
        matplot(time,f_pca[,,3],type="l")
        title(main="f domain: PD 3")
        layout(1)
        cumm_coef = 100*cumsum(s)/sum(s)
        plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")
    }

    return(vfpca)
}
