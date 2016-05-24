if (getRversion() >= "2.15.1") globalVariables(c("x", "y", "f"))
#' @encoding UTF-8
#' @title The Butterfly Curve
#'
#' @description The butterfly curve is a parametric equation discovered by Temple Fay where two functions in a plane produces butterfly-like curves.

#' @param n An integer for background points.
#' @param nb An integer for the butterfly's points.
#' @param title A character vector for plot title.
#' @references
#' Fay, Temple H. (May 1989). The Butterfly Curve. \emph{Amer. Math. Monthly} 96 (5): 442-443. doi:10.2307/2325155.
#' @import ggplot2
#' @export
#'
#' @examples
#' if (interactive()) {
#' butterfly(10, 100, title="10 x 100");
#' butterfly(10, 200, title="10 x 200");
#' butterfly(10, 500, title="100 x 500");
#' butterfly(100, 1000, title="100 x 1000");
#' }
#'
`butterfly` <- function(n=100, nb=500, title=element_blank())
{
s1 <- Sys.time()
t <- seq(0,10*pi,length=nb)

butterfly <- data.frame(x=sin(t)*(exp(1)^cos(t)-2*cos(4*t)-(sin(t/12))^5),
                        y=cos(t)*(exp(1)^cos(t)-2*cos(4*t)-(sin(t/12))^5),
                        s=runif(nb, min=.1, max=10),
                        f=factor(sample(1:10,nb,TRUE)),
                        a=runif(nb,min=.1, max=.4))

points <- data.frame(x=runif(n,-4,4),
                     y=runif(n,-3,5),
                     s=runif(n,min=30, max=50),
                     f=factor(sample(1:10,n,TRUE)),
                     a=runif(n,min=.05, max=.15))

.data <- rbind(butterfly, points)

gplot <- ggplot(.data, aes(x, y, colour=f))
gplot <- gplot + geom_point(alpha=.data$a,size=.data$s)
gplot <- gplot + theme(legend.position="none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.ticks=element_blank(),
          axis.title=element_blank(),
          axis.text =element_blank())
gplot <- gplot + ggtitle(title)

s2 <- Sys.time()
timediff <- c( s2 - s1 )
cat("\n", "Date of Analysis: ",format(Sys.time(), "%a %b %d %Y"), "\n", "Computation time: ",timediff,sep="","\n")
cat("-----------------------------------\n")
return(gplot)
}
