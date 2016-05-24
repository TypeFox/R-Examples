plot.Mort2Dsmooth <-
function(x, 
                              type=c("logrates", "deaths"),
                              palette=c("rainbow",
                                "heat.colors",
                                "terrain.colors",
                                "topo.colors",
                                "cm.colors"), ...){
  object <- x
  type <- match.arg(type)
  palette <- match.arg(palette)
  x <- object$x
  y <- object$y
  Z <- object$Z
  list. <- list(x=x, y=y, type=c("Actual", "Fitted"))
  grid. <- expand.grid(list.)
  Plot <- switch(type,
                 logrates = 1,
                 deaths = 2)
  if(Plot==1){
    ETA <- log(Z) - object$offset
    ETA[object$W==0] <- NA
    ETA.hat <- matrix(MortSmooth_BcoefB(object$Bx,
                                        object$By,
                                        object$coef),
                      length(x),length(y),
                      dimnames = list(x,y))
    grid.$Z <- c(c(ETA), c(ETA.hat))
    my.breaks <- quantile(grid.$Z,
                          prob=seq(0,1,0.1), na.rm=TRUE)
    n.col <- length(my.breaks)-1
    my.col  <- get(palette)(n.col)
  }
  if(Plot==2){
    Z[object$W==0] <- NA
    Z.hat <- object$fitted.values
    grid.$Z <- c(c(Z), c(Z.hat))
    my.breaks <- quantile(grid.$Z,
                          prob=seq(0,1,0.1), na.rm=TRUE)
    n.col <- length(my.breaks)-1
    my.col  <- get(palette)(n.col)
  }
  levelplot(Z ~ y * x | type, grid.,
            layout=c(2,1),
            at=my.breaks, col.regions=my.col, 
            colorkey=list(col=my.col))
}
