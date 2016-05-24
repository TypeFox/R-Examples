
cov.svd <- function(V, chr){
  
  for(i in 0:6){
    a <- 10^(-i)
    svd.obj <- try(svd(V * a), silent = TRUE)
    if(error.try(svd.obj)){
      if(i == 6){
        msg <- paste("SVD error in group", chr)
        stop(msg)
      }else{
        next
      }
    }else{
      svd.obj$d <- svd.obj$d / a
      break
    }
  }
  
  v.svd <- svd.obj$v
  d.svd <- sqrt(abs(svd.obj$d))
  
  rm(svd.obj)
  gc()
  
  U <- d.svd * t(v.svd)
  colnames(U) <- colnames(V)
  
  U
  
}

# 
# cov.svd <- function(V, chr){
#   
#   svd.obj <- try(svd(V), silent = TRUE)
#   if(error.try(svd.obj)){
#     msg <- paste("SVD error in group", chr)
#     stop(msg)
#   }
#   
#   v.svd <- svd.obj$v
#   d.svd <- sqrt(abs(svd.obj$d))
#   if(length(d.svd) > 1){
#     d.svd <- diag(d.svd)
#   }else{
#     dim(d.svd) <- c(1, 1)
#   }
#   rm(svd.obj)
#   gc()
#   
#   U <- d.svd %*% t(v.svd)
#   colnames(U) <- colnames(V)
#   
#   U
#   
# }



