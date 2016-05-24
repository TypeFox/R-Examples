#' metaheurcompare
#'
#' Compare parameter settings
#'
#' @param gridclassobject (GridClass) created by setgrid function in preprocomb package, defaults to examplegrid
#' @param runs (integer) number of runs, defaults to two
#' @param iterations (integer) number of iterations done for a restart
#' @param taboolistlength (integer vector) number of previous solution that can not be revisited, > 1
#' @param initialtemperature (numeric vector) initial propability for acccepting an inferior candidate, between 0 and 1
#' @param tempconst (numeric vector) multiplier for decreasing temperature on each iteration
#' @param reheat (numeric vector) propability of increasing temperature on each iteration
#' @param late (integer vector) location of previous best solution a candidate is compared to, default to 0 for last
#' @examples ## result <- metaheurcompare()
#' @export


metaheurcompare <- function(gridclassobject=examplegrid, runs=2, iterations=10, initialtemperature=c(0.01,1), tempconst=c(1, 0.85),
                            reheat=c(0.01, 0.1), taboolistlength=c(1,2), late=c(0,1)){

  output <- list()
  bestofrun <- numeric()
  histofrun <- list()

  for (u in 1:runs)
  {
    a <- metaheur(gridclassobject, iterations=iterations, initialtemperature=initialtemperature[u],
                  tempconst=tempconst[u], taboolistlength=taboolistlength[u], reheat=reheat[u], late=late[u], stopvalue=0.99, verbose=FALSE)
    b <- a[[3]]
    d <- max(unlist(lapply(b, function(x) x$accuracy)))
    bestofrun[u] <- d
    histofrun[[u]] <- unlist(a[[2]])
  }

  max.length <- max(sapply(histofrun, length))
  l <- lapply(histofrun, function(v) { c(v, rep(NA, max.length-length(v)))})
  res <- t(do.call(rbind, l))
  row.names(res) <- seq(1,nrow(res),1)
  return(res)

}

#' plotsearchpath
#'
#' plot the searh path of metaheurcompare()
#'
#' @param x output of metaheurcompare function
#' @examples ## result <- metaheurcompare()
#' ## plotsearchpath(result)
#' @export

plotsearchpath <- function(x){
  data <- reshape2::melt(x)
  g1 <- ggplot(data, aes(x=Var1, y=value)) + geom_line() + facet_grid(Var2 ~ .) + theme_bw()
  g1
}

#' plotdensity
#'
#' plot the density distribution of classification accuracies of metaheurcompare()
#'
#' @param x output of metaheurcompare function
#' @examples ## result <- metaheurcompare()
#' ## plotdensity(result)
#' @export

plotdensity <- function(x){
  data <- reshape2::melt(x)
  g2 <- ggplot(data, aes(x=value)) + geom_density() + facet_grid(Var2 ~ .) + theme_bw()
  g2
}
















