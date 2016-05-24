do_es <- function(data) {
 # require("corpcor")
  # data: data matrix
  out = list()
  m = nrow(data)  ## number of rows
  p = ncol(data)  ## number of columns
  
  #if (is.cate == F) {
  orgn = apply(data, 2, mean)
  orgn = t(t(orgn))
  # if (is.svd == T) {
  data = t(data)
  #  }
  oner = matrix(1, 1, m)
  #     if (is.svd == F) {
  #       cov.data = (1/n) * (t(data - t(oner) %*% t(orgn)) %*% (data - t(oner) %*% t(orgn)))
  #       eig.res = eigen(cov.data)
  #       #   print('cov.data')
  #       #   print(cov.data)
  #       u = eig.res$vectors
  #       d = eig.res$dues
  #       d = t(t(d))
  #     #  out = list()
  #       out$cov = cov.data
  #     } else {
  #    
  
  cen.mat = data - (orgn) %*% (oner)
  # print(t(cen.mat[1:5,1:6]))
  svd.res = svd(cen.mat)
  u = svd.res$u
  d = svd.res$d
  v = svd.res$v
  #  out = list()
  out$v = u
  #out$cen.mat = cen.mat
  #    }
  #out=list()
  out$m = m
  out$orgn = orgn
  out$u = v
  out$d = d
  out
} 
