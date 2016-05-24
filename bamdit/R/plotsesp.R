#' plotsesp()
#' plot the posterior densities for Se and Sp
#'
#' @param m The matrix used for plotting the diagram.
#' @param binwidth.p Histograms binwidth, default is 0.03.
#' @param CI.level Level of the posterior interval default is 0.95.
#' @seealso \code{\link{metadiag}}.
#' @keywords file
#' @examples
#'
#' ## execute analysis
#' \dontrun{
#' data(ep)
#' m1.ep <- metadiag(ep[,1:4])
#'
#' plotsesp(m = m1.ep)
#' }
#'
#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @import R2jags
#' @import rjags
#' @importFrom stats cor integrate median pnorm qchisq quantile sd
#'


#'@export
plotsesp <- function(m, binwidth.p = 0.03, CI.level = 0.95)
  {

  if(class(m)!="rjags")stop("You have to provide a valid rjags object as fitted model")

  # Patch for the "note" no-visible-binding-for-global-variable
  se.pool <- m$BUGSoutput$sims.list$se.pool
  sp.pool <- m$BUGSoutput$sims.list$sp.pool
  se.new <- m$BUGSoutput$sims.list$se.new
  sp.new <- m$BUGSoutput$sims.list$sp.new

	lo <- (1 - CI.level)/2
  up <- 1 - lo

  p1 <- ggplot(data.frame(se.pool), aes_string(x = "se.pool")) +
          geom_histogram(colour = "black", fill = "blue", binwidth = binwidth.p) +
          xlim(0,1)+
          xlab("Sensitivity (pooled)")+
          geom_vline(xintercept = median(se.pool)) +
          geom_vline(xintercept = quantile(se.pool, prob = c(lo, up)), linetype = "longdash")

  p2 <- ggplot(data.frame(sp.pool), aes_string(x = "sp.pool")) +
          geom_histogram(colour = "black", fill = "blue", binwidth = binwidth.p) +
          xlim(0,1)+
          xlab("Specificity (pooled)")+
          geom_vline(xintercept = median(sp.pool)) +
          geom_vline(xintercept = quantile(sp.pool, prob = c(lo, up)), linetype = "longdash")


  p3 <- ggplot(data.frame(se.new), aes_string(x = "se.new")) +
          geom_histogram(colour = "black", fill = "red", binwidth = binwidth.p) +
          xlim(0,1)+
          xlab("Sensitivity (predictive)")+
          geom_vline(xintercept = median(se.new)) +
          geom_vline(xintercept = quantile(se.new, prob = c(lo, up)), linetype = "longdash")

  p4 <- ggplot(data.frame(sp.new), aes_string(x = "sp.new")) +
          geom_histogram(colour = "black", fill = "red", binwidth = binwidth.p) +
          xlim(0,1)+
          xlab("Specificity (predictive)")+
          geom_vline(xintercept = median(sp.new)) +
          geom_vline(xintercept = quantile(sp.new, prob = c(lo, up)), linetype = "longdash")

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2)))

  vplayout <- function(x, y)
    viewport(layout.pos.row = x,
             layout.pos.col = y)

  print(p1, vp = vplayout(1, 1))
  print(p2, vp = vplayout(1, 2))
  print(p3, vp = vplayout(2, 1))
  print(p4, vp = vplayout(2, 2))
}

