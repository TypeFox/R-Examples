#' Horizontal Functional Principal Component Analysis
#'
#' This function calculates vertical functional principal component analysis
#' on aligned data
#'
#' @param gam matrix (\eqn{N} x \eqn{N}) of \eqn{M} of warping functions with \eqn{N} time samples
#' @param no number of prinicpal components to extract
#' @param showplot show plots of prinipal directions (default = T)
#' @return Returns a list containing \item{gam_pca}{warping functions principal directions}
#' \item{psi_pca}{srvf principal directions}
#' \item{latent}{latent values}
#' \item{U}{eigenvectors}
#' \item{vec}{shooting vectors}
#' \item{mu}{Karcher Mean}
#' @keywords srvf alignment
#' @references Tucker, J. D., Wu, W., Srivastava, A.,
#'  Generative Models for Function Data using Phase and Amplitude Separation,
#'  Computational Statistics and Data Analysis (2012), 10.1016/j.csda.2012.12.001.
#' @export
#' @examples
#' data("simu_warp")
#' hfpca = horizFPCA(simu_warp$gam,no = 3)
horizFPCA <- function(gam,no,showplot = TRUE){
    tmp = SqrtMean(gam)
    vec = tmp$vec
    mu = tmp$mu
    # Parameters
    tau = 1:5 # -2, -1, 0, 1, 2 std from the mean
    no_pca = 1:no

    # TFPCA
    K = cov(t(vec)) #out$sigma

    out = svd(K)
    s = out$d
    U = out$u
    TT = nrow(vec) + 1
    vm = rowMeans(vec)

    gam_pca = array(0,dim=c(length(tau),length(mu)+1,no))
    psi_pca = array(0,dim=c(length(tau),length(mu),no))
    for (j in no_pca){
        cnt = 1
        for (k in tau){
            v = (k-3)*sqrt(s[j])*U[,j]
            vn = pvecnorm(v,2)/sqrt(TT)
            if (vn < 0.0001){
                psi_pca[k,,j] = mu
            }else{
                psi_pca[k,,j] = cos(vn)*mu+sin(vn)*v/vn
            }
            tmp = rep(0,TT)
            tmp[2:TT] = cumsum(psi_pca[k,,j]*psi_pca[k,,j])
            gam_pca[k,,j] = tmp/TT
            cnt = cnt + 1
        }
    }

    hfpca = list()
    hfpca$gam_pca = gam_pca
    hfpca$psi_pca = psi_pca
    hfpca$latent = s
    hfpca$U = U
    hfpca$vec = vec
    hfpca$mu = mu

    if (showplot){
        layout(matrix(c(1,2,3), 1, 3, byrow = TRUE))
        matplot(seq(0,1,len=TT),t(gam_pca[,,1]),type="l",xlab = "t",ylab = "t")
        title(main="PD 1")
        matplot(seq(0,1,len=TT),t(gam_pca[,,2]),type="l",xlab = "t",ylab = "t")
        title(main="PD 2")
        matplot(seq(0,1,len=TT),t(gam_pca[,,3]),type="l",xlab = "t",ylab = "t")
        title(main="PD 3")
        layout(1)
        cumm_coef = 100*cumsum(s)/sum(s)
        plot(cumm_coef,type="l",col="blue",main="Coefficient Cumulative Percentage", ylab = "Percentage")
    }
    return(hfpca)
}
