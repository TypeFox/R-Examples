#' Plot fitted population pdfs
#' 
#' Plot fitted population pdfs at combinations of two different values of the population mean and variance.
#' 
#' Four pdfs are plotted, using each combination of the \code{alpha}/2 and 1-\code{alpha}/2
#' quantiles of the fitted distributions for the population median and standard deviation
#' 
#' @param medianfit The output of a \code{fitdist} command following elicitation
#'  of the expert's beliefs about the population median.
#' @param precisionfit The output of a \code{fitdist} command following elicitation
#'  of the expert's beliefs about the population precision.
#'  
#' @param alpha Value between 0 and 1 to determine choice of means and variances used in plots
#' @param tails Value between 0 and 1 to determine the tail area shown in the pdf plots
#' @param lower lower limit on the x-axis for plotting.
#' @param upper upper limit on the x-axis for plotting.
#' @param d The fitted distribution for the population median. Can be one of "normal",
#'  "lognormal" or "best", where "best" will select the best fitting out of 
#'  normal and lognormal.
#' @param n.x The number of points on the x-axis at which the pdf is plotted.
#' @param fontsize Font size used in the plots.
#'  
#' @return A plot and a list, containing
#' \item{mu}{The two population mean values used in the plots.}
#' \item{sigma}{The two population standard deviation values used in the plots.}
#'  
#'  
#' @references \code{multiplot} function obtained from
#'  \url{http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/}
#'  
#' @examples 
#' \dontrun{
#' prfit <- fitprecision(interval = c(60, 70), propvals = c(0.2, 0.4), trans = "log")
#' medianfit <- fitdist(vals = c(50, 60, 70), probs = c(0.05, 0.5,  0.95), lower = 0)
#' pdfplots(medianfit, prfit, alpha = 0.01)
#'  }
#' @import ggplot2
#' @export

pdfplots <- function(medianfit, precisionfit, 
                     alpha = 0.05, tails = 0.05, 
                     lower = NA, upper = NA, n.x = 100, 
                     d = "best",
                     fontsize = 18){
  
  mediandist <- getmediandist(medianfit, d)
  
  a<-precisionfit$Gamma[[1]]
  b<-precisionfit$Gamma[[2]]
  
  quantiles <- c(alpha / 2, 1 - alpha /2)
  
  f <- getdists(precisionfit$transform)
  
  lim <- getlimits(lower, upper, f, mediandist, precisionfit)
  
  mu <- f$trans(mediandist$quan(quantiles, mediandist$m ,mediandist$s))
  sigma <- sort(sqrt(1 / qgamma(quantiles, a, b)))
  
  X <- seq(from = lim$lower, to = lim$upper, length = n.x)
  
  f11 <- f12 <- f21 <- f22 <- NULL # hack to avoid R CMD check NOTE 
  
  df <- data.frame(X, f11 = f$dens(X, mu[1], sigma[1]), 
                   f12 = f$dens(X, mu[1], sigma[2]),
                   f21 = f$dens(X, mu[2], sigma[1]), 
                   f22 = f$dens(X, mu[2], sigma[2]))
  
  dmax <- max(df[,2:5])
  
  theme_set(theme_grey(base_size = fontsize))
  
  pcore <- ggplot(df, aes(x=X, y=f11)) + expand_limits(y = c(0, dmax)) +
    labs(y = "")
  p1 <-  pcore + geom_line() + labs (title = bquote(list(mu[.(alpha/2)], ~ sigma[.(alpha/2)])))
  p2 <- pcore + geom_line(aes(x=X, y=f12)) + labs (title = bquote(list(mu[.(alpha/2)], ~ sigma[.(1 - alpha/2)])))
  p3 <- pcore + geom_line(aes(x=X, y=f21)) + labs (title = bquote(list(mu[.(1 - alpha/2)], ~ sigma[.(alpha/2)])))
  p4 <- pcore + geom_line(aes(x=X, y=f22)) + labs (title = bquote(list(mu[.(1 - alpha/2)], ~ sigma[.(1 - alpha/2)])))
  
  if(!is.na(tails)){
    
    xl <- dl <- xu <- du <- NULL # hack to avoid R CMD check NOTE
    
    df <- taildensities(mu[1], sigma[1], tails, n.x, lim$lower, lim$upper, f$dens, f$quan)
    p1 <- p1 + geom_area(data = df, aes(x=xl, y = dl), fill="red", alpha=0.5) +
      geom_area(data = df, aes(x=xu, y = du), fill="red", alpha=0.5)
    
    df <- taildensities(mu[1], sigma[2], tails, n.x, lim$lower, lim$upper, f$dens, f$quan)
    p2 <- p2 + geom_area(data = df, aes(x=xl, y = dl), fill="red", alpha=0.5) +
      geom_area(data = df, aes(x=xu, y = du), fill="red", alpha=0.5)
    
    df <- taildensities(mu[2], sigma[1], tails, n.x, lim$lower, lim$upper, f$dens, f$quan)
    p3 <- p3 + geom_area(data = df, aes(x=xl, y = dl), fill="red", alpha=0.5) +
      geom_area(data = df, aes(x=xu, y = du), fill="red", alpha=0.5)
    
    df <- taildensities(mu[2], sigma[2], tails, n.x, lim$lower, lim$upper, f$dens, f$quan)
    p4 <- p4 + geom_area(data = df, aes(x=xl, y = dl), fill="red", alpha=0.5) +
      geom_area(data = df, aes(x=xu, y = du), fill="red", alpha=0.5)
    
  }
  multiplot(p1, p2, p3, p4, cols = 2)
  
  list(mu = mu, sigma = sigma)
}