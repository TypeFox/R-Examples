# This is R source code for function 'threeStageParSel', in the
# R package "image".
# Date: Aug 30, 2015
# Creator: Yicheng Kang

threeStageParSel = function(image, bandwidth, edge1, edge2, nboot, blur=FALSE){
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
  edge1 = matrix(as.integer(edge1), ncol=n1)
  edge2 = matrix(as.integer(edge2), ncol=n1)
  n_band = length(bandwidth)
  if (blur == FALSE) {
      out = .Fortran('denoise_3stage_bandwidth', n=as.integer(n1-1), obsImg=z, nband=n_band,
          bandwidth=as.integer(bandwidth), edge1=edge1, edge2=edge2, cv=rep(as.double(0), n_band))
      k.cv = out$cv
      band_sel = mean(bandwidth[k.cv ==  min(k.cv)])
      out.mat = matrix(0, ncol=n_band, nrow=1)
      out.mat[1, ] = k.cv
      colnames(out.mat) = paste('bandwidth=', bandwidth, sep='')
      rownames(out.mat) = "CV-score"

  }
  else {
      if(missing(nboot))
          stop("nboot must be specified when blur is TRUE.")
      if (length(nboot) > 1)
          stop("nboot must be an integer number.")
      n_boot = as.integer(nboot)
      out = .Fortran('deblur_3stage_bandwidth', n=as.integer(n1-1), obsImg=z, nband=n_band,
          bandwidth=as.integer(bandwidth), edge1=edge1, edge2=edge2, nboot=n_boot,
          msecv=rep(as.double(0), n_band))
      k.msecv = out$msecv
      band_sel = mean(bandwidth[k.msecv ==  min(k.msecv)])
      out.mat = matrix(0, ncol=n_band, nrow=1)
      out.mat[1, ] = k.msecv
      colnames(out.mat) = paste('bandwidth=', bandwidth, sep='')
      rownames(out.mat) = "MSECV-score"
  }
  
  print(paste('The selected bandwidth is', band_sel))
  return(list(output_matrix=out.mat, selected_bandwidth=band_sel))
}

