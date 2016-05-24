# This is R source code for function 'surfaceCluster', in the
# R package "image".
# Date: May 3, 2013
# Creator: Yicheng Kang

surfaceCluster = function(image, bandwidth, sig.level,
  blur = FALSE, plot = FALSE){
  if (!is.matrix(image))
    stop("image data must be a matrix")
  else n1 = dim(image)[1]
       n2 = dim(image)[2]
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
  if (length(bandwidth) > 1) {
    stop('bandwidth must be a positive integer')
  }
  if (n1 + 2 * bandwidth + 2 > 600)
    stop("some choice of bandwidth or the resolution of the
image is too large")
  if (!is.numeric(sig.level) | abs(sig.level - 0.5) > 0.5)
    { stop('sig.level must be a number between 0 and 1') }
  n1 = dim(image)[1]
  z = matrix(as.double(image), ncol = n1)
  zq = as.double(qnorm(sig.level))
  jp.llk = JPLLK_surface(z, 2:7)
  fitted = jp.llk$fitted
  resid = jp.llk$resid
  mmt2 = (jp.llk$sigma)^2
  mmt4 = mean(resid^4)
  k = as.integer(bandwidth)
  if (blur == FALSE) {
    out = .Fortran('cluster_denoise', n = as.integer(n1 - 1),
      obsImg = z, bandwidth = k, mmt2 = as.double(mmt2), mmt4 =
      as.double(mmt4), zq = zq, estImg = z)
  }
  else {
    out = .Fortran('cluster_deblur', n = as.integer(n1 - 1),
            obsImg = z, bandwidth = k, mmt2 = as.double(mmt2), mmt4 =
            as.double(mmt4), zq = zq, estImg = z)
  }
  if (plot == FALSE) {
    return(list(estImg = out$estImg, sigma2 = out$mmt2,
                mmt4 = out$mmt4))
  }
  else { x = seq(0, 1, length = n1); y = x
         image(x, y, out$estImg, zlim = c(0, 255), col = gray((0:255)/
                                                     255))
         return(list(estImg = out$estImg, sigma2 = out$mmt2,
                     mmt4 = out$mmt4))
       }
}
