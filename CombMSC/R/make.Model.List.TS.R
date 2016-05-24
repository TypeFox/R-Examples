`make.Model.List.TS` <-
function(max.Orders, period = 12){
  make.TSM <- function(vec){
    obj <- list()
    obj$model <- list(AR = vec[1], "I" = vec[2], MA = vec[3])
    obj$seasonal <- list()
    obj$seasonal$model <- list(AR = vec[4], "I" = vec[5], MA = vec[6])
    obj$seasonal$period = period
    class(obj) <- "tsm"
    obj
  }
  if(!is.numeric(max.Orders) || length(max.Orders) != 6) stop("max.Orders must be a numeric vector of length 6.")
  mat <- list()
  mat[[1]] <- rep(0:max.Orders[1], times = 1, each = prod(max.Orders[2:6]+1))
  for(j in 2:5){
    mat[[j]] <- rep(0:max.Orders[j], each = prod(max.Orders[(j+1):6]+1), times = prod(max.Orders[1:(j-1)]+1))
  }
  mat[[6]] <- rep(0:max.Orders[6], times = prod(max.Orders[1:5]+1), each = 1)
  mat <- unlist(mat)
  dim(mat) <- c(prod(max.Orders+1), 6)
  final <- apply(mat, 1, make.TSM)
  final
}

