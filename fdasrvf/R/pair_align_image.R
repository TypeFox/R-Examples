#' Pairwise align two images

#' This function aligns to images using the q-map framework
#'
#' @param I1 reference image
#' @param I2 image to warp
#' @param M number of basis elements (default=5)
#' @param ortho orthonormalize basis (default=TRUE)
#' @param basis_type ("t","s","i","o"; default="t")
#' @param resizei resize image (default=TRUE)
#' @param N size of resized image (default=64)
#' @param stepsize gradient stepsize (default=1e-5)
#' @param itermax maximum number of iterations (default=1000)
#' @return Returns a list containing \item{Inew}{aligned I2}
#' \item{gam}{warping function}
#' @keywords image alignment
#' @references Q. Xie, S. Kurtek, E. Klassen, G. E. Christensen and A. Srivastava. Metric-based pairwise and multiple image registration. IEEE European Conference on Computer Vision (ECCV), September, 2014
#' @export
pair_align_image <- function(I1, I2, M=5, ortho=TRUE, basis_type="t", resizei=TRUE, N=64, stepsize=1e-5, itermax=1000){
    m = dim(I1)[1]
    n = dim(I1)[2]
    F1 = array(0,dim=c(m,n,2))
    m1 = dim(I2)[1]
    n1 = dim(I2)[1]
    F2 = array(0,dim=c(m1,n1,2))

    # Take Gradient-------------------------------------------------------------
    out = gradient2(I1,1./(m-1), 1./(n-1))
    F1[,,1] = out$dxdu
    F1[,,2] = out$dydv
    out = gradient2(I2,1./(m1-1), 1./(n1-1))
    F2[,,1] = out$dxdu
    F2[,,2] = out$dydv

    # Resize Data and Center----------------------------------------------------
    if (resizei){
        if ((N>m) || (N>n)){
            cat("Not resizing, N is larger than image size")
        } else {
            xlim = c(1,m)
            ylim = c(1,n)
            dx = (m-1)/(N-1)
            dy = (n-1)/(N-1)
            F1a = array(0,dim=c(N,N,2))
            if (requireNamespace("akima", quietly = TRUE)) {
                F1a[,,1] = akima::bicubic.grid(1:m,1:n,F1[,,1],xlim,ylim,dx,dy)$z
                F1a[,,2] = akima::bicubic.grid(1:m,1:n,F1[,,2],xlim,ylim,dx,dy)$z
            } else {
                grid.list<- list(x=seq(1,m,length.out=N), y=seq(1,n,length.out=N))
                obj<-list(x=1:m, y=1:n, z=F1[,,1])
                F1a[,,1] = interp.surface.grid(obj, grid.list)$z
                obj<-list(x=1:m, y=1:n, z=F1[,,2])
                F1a[,,2] = interp.surface.grid(obj, grid.list)$z
            }
            F1 = F1a

            xlim = c(1,m1)
            ylim = c(1,n1)
            dx = (m1-1)/(N-1)
            dy = (n1-1)/(N-1)
            F2a = array(0,dim=c(N,N,2))
            if (requireNamespace("akima", quietly = TRUE)) {
              F2a[,,1] = akima::bicubic.grid(1:m1,1:n1,F2[,,1],xlim,ylim,dx,dy)$z
              F2a[,,2] = akima::bicubic.grid(1:m1,1:n1,F2[,,2],xlim,ylim,dx,dy)$z
            } else {
              grid.list<- list(x=seq(1,m1,length.out=N), y=seq(1,n1,length.out=N))
              obj<-list(x=1:m1, y=1:n1, z=F2[,,1])
              F2a[,,1] = interp.surface.grid(obj, grid.list)$z
              obj<-list(x=1:m1, y=1:n1, z=F2[,,2])
              F2a[,,2] = interp.surface.grid(obj, grid.list)$z
            }
            F2 = F2a
        }
    }

    F1 = F1 - min(F1)
    F1 = F1 / max(F1)
    F2 = F2 - min(F2)
    F2 = F2 / max(F2)

    # Gernerate basis-----------------------------------------------------------
    out = run_basis(F1, M, basis_type, ortho)
    b = out$b
    gamid = out$gamid
    gamp = gamid

    out = reparam_image(F1, F2, gamp, b, stepsize=stepsize, itermax=itermax)

    I2_new = apply_gam_to_imag(I2,out$gam)

    return(list(I2_new=I2_new, gam=out$gam))
}
