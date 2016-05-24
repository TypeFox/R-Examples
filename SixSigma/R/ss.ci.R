if(getRversion() >= '2.15.1') utils::globalVariables(c("..density.."))
#' Confidence Interval for the mean
#' 
#' Computes a confidence interval for the mean of the variable (parameter
#' or feature of the process), and prints the data, a histogram with a density line,
#' the result of the Shapiro-Wilks normality test and a quantile-quantile plot.
#' 
#' When the population variance is known, or the size is greater than 30,
#' it uses z statistic. Otherwise, it is uses t statistic.\cr
#' If the sample size is lower than 30, a warning is displayed so as to
#' verify normality.
#' 
#' @param x         A numeric vector with the variable data
#' @param sigma2    The population variance, if known
#' @param alpha     The eqn{\\alpha} error used to compute the \eqn{100*(1-\\alpha)\%} confidence interval
#' @param data      The data frame containing the vector
#' @param xname     The name of the variable to be shown in the graph
#' @param approx.z  If TRUE it uses z statistic instead of t when sigma is unknown and sample size 
#' is greater than 30. The default is FALSE, change only if you want to compare with
#' results obtained with the old-fashioned method mentioned in some books.
#' @param main      The main title for the graph
#' @param digits    Significant digits for output
#' @param sub       The subtitle for the graph (recommended: six sigma project name)
#' @param ss.col    A vector with colors
#' @return 
#'   The confidence Interval.\cr
#'   A graph with the figures, the Shapiro-Wilks test,  and a histogram.
#' 
#' @references 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#' 
#' @author EL Cano
#' 
#' @note 
#' Thanks to the kind comments and suggestions from the anonymous reviewer 
#' of a tentative article.
#' 
#' @seealso 
#' \code{\link{ss.data.rr}}
#' @examples
#' ss.ci(len, data=ss.data.strings, alpha = 0.05,
#'   sub = "Guitar Strings Test | String Length", 
#'   xname = "Length")
#' @keywords confidence interval normality test mean
#' @export
ss.ci<-function(x, sigma2 = NA, alpha = 0.05, data = NA, 
                xname = "x", approx.z = FALSE, main = "Confidence Interval for the Mean", 
                digits = 3,
                sub = "", ss.col = c("#666666", "#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")){
  if (is.data.frame(data)){
    x <- data[[deparse(substitute(x))]]
  }
  na <- length(which(is.na(x)))
  if (na > 0) { cat(na, " missing values were ommitted\n")}
  m <- mean(x, na.rm = TRUE)
  n <- length(x) - na
  s <- ifelse(is.numeric(sigma2), 
              sqrt(sigma2),
              sd(x, na.rm = TRUE))
  if (is.numeric(sigma2) | approx.z == TRUE){
    st <- qnorm(1 - (alpha/2))
    st.dist <- c("z")
  }
  else{
    if (n < 30) {
      warning("\nThe sample size is lower than 30. Check Normality\n\n")
    }
    st <- qt(1-(alpha/2), n-1)
    st.dist <- c("t")
  }
  dist <- st * (s/sqrt(n))
  
  cat("\tMean = ", round(m, digits), "; sd = ", round(s, digits), "\n", sep = "")
  cat("\t", (1-alpha)*100, "% Confidence Interval= ", 
      round(m-dist, digits), " to ", round(m+dist, digits),"\n\n", sep = "")
  ci <- c(m-dist, m+dist)
  names(ci) <- c("LL", "UL")
  ##Canvas-container
  .ss.prepCanvas(main,sub, ss.col)
  
  ##figures
  vp.figures <- grid::viewport(name = "figures", 
                               x = 0, 
                               width = 1, 
                               height = unit(8, "lines"),
                               y = 1, 
                               just = c("left", "top"),
                               layout = grid::grid.layout(1, 2, widths = c(0.4, 0.6))						 )
  grid::pushViewport(vp.figures)
  
  vp.figures1 <- grid::viewport(name = "figures1", 
                                layout.pos.col = 1, 
                                layout.pos.row = 1)
  grid::pushViewport(vp.figures1)
  grid:: grid.roundrect(height = unit(7, "lines"), 
                        width = 0.95,
                        gp = grid::gpar(fill = ss.col[5], col = ss.col[2], lwd=2))
  grid:: grid.text("Mean:\nStdDev:\nn:\nMissing:", just = "left",
                   x = unit(1, "npc") - unit(5.5, "cm"), 
                   gp = grid::gpar(fontface = c("bold")))
  grid:: grid.text(paste(round(m, digits), "\n", round(s, digits), "\n", n,
                         "\n", na, sep = ""), just = "right", 
                   x = unit(1, "npc") - unit(1, "cm"))
  
  
  grid::popViewport()
  vp.figures2 <- grid::viewport(name = "figures2", layout.pos.col = 2,
                                layout.pos.row = 1)
  grid::pushViewport(vp.figures2)
  grid:: grid.roundrect(height = unit(7, "lines"), 
                        width = 0.95,
                        gp = grid::gpar(fill = ss.col[5], col = ss.col[2], lwd = 2))
  
  grid:: grid.text(paste((1-alpha)*100, "% CI:\nP-Var:\n",
                         st.dist, ":", sep = ""),
                   just = "left",
                   x = unit(0,"npc") + unit(1,"cm"), 
                   gp = grid::gpar(fontface=c("bold")))
  grid:: grid.text(paste("[", round(ci[1], digits), ", ", round(ci[2], digits),
                         "]\n", ifelse(is.numeric(sigma2), sigma2, "unknown"),
                         "\n", round(st, digits),
                         sep = ""), just = "right", 
                   x = unit(0,"npc") + unit(7.5,"cm"))
  grid::popViewport()
  grid::popViewport()
  
  #graph
  vp.graph <- grid::viewport(name = "graph", 
                             y = 0, 
                             width = 0.95,
                             height = unit(1, "npc") - unit(8, "lines"),
                             just = c("center", "bottom"),
                             layout = grid::grid.layout(1, 2, 
                                                        widths = unit(c(1, 6), c("null", "cm"))))
  grid::pushViewport(vp.graph)
  
  vp.test <- grid::viewport(name = "test", layout.pos.row = 1, layout.pos.col = 2)
  grid::pushViewport(vp.test)
  grid::grid.rect()
  grid:: grid.roundrect(height = unit(6, "lines"),
                        width = 0.9, 
                        y = unit(1, "npc") + unit(-1, "lines"),
                        just = "top",
                        gp = grid::gpar(fill = ss.col[5], col = ss.col[2], lwd = 2))
  
  grid:: grid.text("Shapiro-Wilks\nNormality Test\n", 
                   y = unit(1, "npc") - unit(3, "lines"),
                   gp = grid::gpar(fontface = c("bold")))
  pval <- shapiro.test(x)[2]$p.value 
  if (pval < 0.05){
    warning("Sample data is non-normal.")
  }
  grid:: grid.text(paste(round(shapiro.test(x)[1]$statistic, digits), "\n"),
                   y = unit(1, "npc") - unit(5, "lines"))
  grid:: grid.text(paste("p-value:", round(pval, digits), "\n"),
                   y = unit(1, "npc") - unit(6, "lines"))
  
  vp.qq <- grid::viewport(name="qqp", 
                          x = 0.5, y=0.25,
                          height = unit(0.6,"npc"))
  grid::pushViewport(vp.qq)
  qqp <- qplot(sample = x
               # , stat="qq"
  ) + 
    xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank()) +
    ggtitle("Normal q-q Plot")
  
  print(qqp,newpage=FALSE)  
  grid::popViewport()
  grid::popViewport()
  
  vp.hist <- grid::viewport(name = "hist", 
                            layout.pos.row = 1,
                            layout.pos.col = 1)
  grid::pushViewport(vp.hist)
  grid::grid.rect()
  
  ggdata <- reshape2::melt(x)
  nbins <- nclass.Sturges(x)
  qqp <- ggplot(ggdata, aes(x = value))
  myhist <- qqp + 
    geom_histogram(aes(y = ..density..), 
                   bins = nbins,
                   # binwidth = binw,
                   fill = "white",
                   col = "gray"
                   # , stat = "bin"
    ) +
    xlab(paste("Value of", deparse(substitute(x)))) +
    ggtitle("Histogram & Density Plot") +
    stat_density(geom = "path", 
                 position = "identity", 
                 # binwidth = binw,
                 size = 1) 
  suppressWarnings(
    print(myhist, newpage=FALSE)
  )
  
  return (ci)
}
