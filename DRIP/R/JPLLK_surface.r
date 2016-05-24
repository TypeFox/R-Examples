# This is R source code for function 'JPLLK_surface', in the
# R package "image".
# Date: April 25, 2013
# Creator: Yicheng Kang

JPLLK_surface = function(image, bandwidth, plot = FALSE){
  if (!is.matrix(image))
    stop("image data must be a matrix")
  else n1 = dim(image)[1]
       n2 = dim(image)[2]
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
  if (n1 + 2 * max(bandwidth) + 2 > 600)
    stop("some choice of bandwidth or the resolution of the
image is too large")
  n1 = dim(image)[1]
  z = matrix(as.double(image), ncol = n1)
  n_band = length(bandwidth)
  out = .Fortran('JP_LLK_CV', n=as.integer(n1-1), obsImg=z, nband=n_band,
      bandwidth=as.integer(bandwidth), cv=rep(as.double(0), n_band))
  k.cv = out$cv
  cv.band = mean(bandwidth[k.cv ==  min(k.cv)])
  jp.llk = .Fortran('JP_LLK_fit', n = as.integer(n1 - 1),
    obsImg = z, bandwidth = as.integer(cv.band), fitted = z, resid
    = z, sigma = as.double(0))
  if (plot == FALSE) {
    return(list(fitted = jp.llk$fitted, cv.band = cv.band, resid
                = jp.llk$resid, sigma = jp.llk$sigma))
  }
  else { x = seq(0, 1, length = n1); y = x
         image(x, y, jp.llk$fitted, zlim = c(0, 255), col =
               gray((0:255)/255))
         return(list(fitted = jp.llk$fitted, cv.band = cv.band, resid
                = jp.llk$resid, sigma = jp.llk$sigma))
  }
}

