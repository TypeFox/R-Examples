#' Fit Mallows model N times and select most likely model.
#' The EM algorithm to fit Multi-Modal Mallows' models is prone to getting
#' stuck in local maxima, so we run it several times and selec the best one.
#'
#' @param N number of times to run the model
#' @param iter maximum number of iterations for each run
#' @param G Number of cluster centers
#' @param datas data set to fit
#' @return best fitting model.
BestFit <- function(datas, N, iter, G) {
  mod <- list()
  for (i in 1:N) {
    mod[[i]] <- Mallows(datas, iter = iter, G = G)
  }
  best.likes <- unlist(lapply(mod, function(i) i[[5]][length(which(i[[5]] != 0))]))
  best <- which.max(best.likes)
  out <- mod[[best]]
  return(out)
}