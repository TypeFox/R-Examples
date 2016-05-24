#' Plot an \code{mlmc.test} object
#'
#' Produces diagnostic plots on the result of an \code{\link{mlmc.test}} function call.
#'
#' @param x an \code{mlmc.test} object as produced by a call to the \code{\link{mlmc.test}} function.
#' @param which a vector of strings specifying which plots to produce, or \code{"all"} to do all diagnostic plots.  The options are: \describe{
#'   \item{\code{"var"} = \eqn{log_2} of variance against level;}{}
#'   \item{\code{"mean"} = \eqn{log_2} of mean against level;}{}
#'   \item{\code{"consis"} = consistency against level;}{}
#'   \item{\code{"kurt"} = kurtosis against level;}{}
#'   \item{\code{"Nl"} = \eqn{log_2} of number of samples against level;}{}
#'   \item{\code{"cost"} = \eqn{log_10} of cost against \eqn{log_10} of epsilon (accuracy).}{}
#' }
#' @param cols the number of columns across to plot to override the default value.
#' @param ... additional arguments which are passed on to plotting functions.
#'
#' @author Louis Aslett <aslett@stats.ox.ac.uk>
#'
#' @examples
#' \dontrun{
#' tst <- mlmc.test(opre_l, M=4, N=2000000,
#'                  L=5, N0=1000,
#'                  eps.v=c(0.005, 0.01, 0.02, 0.05, 0.1),
#'                  Lmin=2, Lmax=6, option=1)
#' tst
#' plot(tst)
#' }
#'
#' @importFrom ggplot2 ggplot aes_string geom_point geom_line xlab ylab scale_x_log10 scale_y_log10 annotation_logticks
#' @export
plot.mlmc.test <- function(x, which="all", cols=NA, ...) {
  if(length(which)==1 && which=="all") {
    which <- c("var", "mean", "consis", "kurt", "Nl", "cost")
  }
  p <- list()
  if("var" %in% which) {
    p <- c(p, list(
      ggplot(data.frame(l=rep(0:x$L, 2),
                        var=c(log2(x$var1), log2(x$var2)),
                        Method=c(rep("MLMC", x$L+1), rep("MC", x$L+1)))) +
        geom_point(aes_string(x="l", y="var", colour="Method")) +
        geom_line(aes_string(x="l", y="var", colour="Method", linetype="Method")) +
        xlab("Level") +
        ylab(expression(log[2](Variance)))
    ))
  }
  if("mean" %in% which) {
    p <- c(p, list(
      ggplot(data.frame(l=rep(0:x$L, 2),
                        mean=c(log2(x$del1), log2(x$del2)),
                        Method=c(rep("MLMC", x$L+1), rep("MC", x$L+1)))) +
        geom_point(aes_string(x="l", y="mean", colour="Method")) +
        geom_line(aes_string(x="l", y="mean", colour="Method", linetype="Method")) +
        xlab("Level") +
        ylab(expression(log[2](Mean)))
    ))
  }
  if("consis" %in% which) {
    p <- c(p, list(
      ggplot(data.frame(l=0:x$L,
                        consis=x$chk1)) +
        geom_point(aes_string(x="l", y="consis")) +
        geom_line(aes_string(x="l", y="consis")) +
        xlab("Level") +
        ylab("Consistency check")
    ))
  }
  if("kurt" %in% which) {
    p <- c(p, list(
      ggplot(data.frame(l=0:x$L,
                        kurt=x$kur1)) +
        geom_point(aes_string(x="l", y="kurt")) +
        geom_line(aes_string(x="l", y="kurt")) +
        xlab("Level") +
        ylab("Kurtosis")
    ))
  }
  if("Nl" %in% which) {
    p <- c(p, list(
      ggplot(data.frame(l=unlist(lapply(sapply(x$Nl, length), seq)),
                        Nl=log2(unlist(x$Nl)),
                        Epsilon=as.factor(rep(x$eps.v, times=sapply(x$Nl, length))))) +
        geom_point(aes_string(x="l", y="Nl", colour="Epsilon")) +
        geom_line(aes_string(x="l", y="Nl", colour="Epsilon", linetype="Epsilon")) +
        xlab("Level") +
        ylab(expression(log[2](N[l])))
    ))
  }
  if("cost" %in% which) {
    p <- c(p, list(
      ggplot(data.frame(eps=rep(x$eps.v, 2),
                        cost=c(x$eps.v^2*x$mlmc_cost, x$eps.v^2*x$std_cost),
                        Method=c(rep("MLMC", length(x$eps.v)), rep("MC", length(x$eps.v))))) +
        geom_point(aes_string(x="eps", y="cost", colour="Method")) +
        geom_line(aes_string(x="eps", y="cost", colour="Method", linetype="Method")) +
        xlab(expression(log[10](epsilon))) +
        ylab(expression(log[10](epsilon^2*Cost))) +
        scale_x_log10() +
        scale_y_log10() +
        annotation_logticks()
    ))
  }
  if(is.na(cols)) {
    if(length(p) <= 3)
      cols <- length(p)
    if(length(p) == 4)
      cols <- 2
    if(length(p) > 4)
      cols <- 3
  }

  p <- c(p, list(cols=cols))
  do.call(multiplot, p)
}
