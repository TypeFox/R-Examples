## grid.control.R

grid.control <- function(circ.n = 4,
                         circ.r = NULL,
                         circ.col = "gray30",
                         circ.lwd = 0.5,
                         circ.cex = 0.8,
                         circ.between = 0.3,
                         circ.dir = pi/16*9,
                         circ.sub.n = NULL,
                         circ.sub.r = NULL,
                         circ.sub.col = "gray70",
                         circ.sub.lwd = 0.5,
                         cyclVar.lab =
                         c("N","NE","E","SE","S","SW","W","NW"),
                         cyclVar.cex = 1.2,
                         cyclVar.between = 0,
                         cyclVar.centered = TRUE,
                         ray.lim = NULL,
                         ray.n = 8)
    ## helper function for plot.rose
    ## Author: Rene Locher
    ## Version: 2007-10-02

  {
    return(list(ray =
                list(n = ray.n,
                     lim = ray.lim),
                circ =
                list(n = circ.n,
                     r = circ.r,
                     value = NULL,
                     col = circ.col,
                     lwd = circ.lwd,
                     cex = circ.cex,
                     dir = circ.dir,
                     between = circ.between,
                     sub = list(
                       plot=FALSE,
                       n = circ.sub.n,
                       r = circ.sub.r,
                       col = circ.sub.col,
                       lwd = circ.sub.lwd)
                     ),
                cyclVar =
                list(lab = cyclVar.lab,
                     cex = cyclVar.cex,
                     between = cyclVar.between,
                     centered = cyclVar.centered)))
  } ## grid.control


