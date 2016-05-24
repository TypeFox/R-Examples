##############################################################################
#
# calculate composite log likelihood for a potts model.
#
##############################################################################
composite.ll <- function(theta, t_stat, t_cache=NULL, fapply=lapply) {
  tot <- 0
  if (is.null(t_cache)) {
    stop("!cache_t not implemented yet in composite.ll!")
  } else {
    tot_list <- fapply(t_cache, function(arr) {
                                        # subtact base case
      tmp <- t(apply(arr, 1, function(r) r - arr[1,]))
      (t_stat[-1] - arr[1,]) %*% theta - log(sum(exp(tmp %*% theta)))
    })
    sum(unlist(tot_list))
  }
}

##############################################################################
#
# calculate gradient of the composite log likelihood for a potts
# model.
#
##############################################################################
gr.composite.ll <- function(theta, t_stat, t_cache=NULL, fapply=lapply) {
  tot <- rep(0,length(theta))
  if (is.null(t_cache)) {
    stop("!cache_t not implemented yet in composite.ll!")
  } else {
    tot_list <- fapply(t_cache, function(arr) {
      tmp <- t( t(arr) - arr[1,] )
      rtmp <- as.vector( exp( tmp %*% theta ) )
      t_stat[-1] - arr[1,] - colSums( tmp * rtmp / sum(rtmp) )
    })
    colSums(matrix(unlist(tot_list), ncol=length(theta), byrow=TRUE))
  }
}
