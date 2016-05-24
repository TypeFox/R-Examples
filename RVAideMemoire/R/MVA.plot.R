MVA.plot <- function(x,type=c("scores","loadings","correlations","biplot","pairs",
  "trajectories"),...) {
  type <- match.arg(type)
  x <- MVA.ident(x)
  f <- switch(type,scores=MVA.scoreplot,loadings=MVA.loadplot,correlations=MVA.corplot,
    biplot=MVA.biplot,pairs=MVA.pairplot,trajectories=MVA.trajplot)
  f(x,...)
}


