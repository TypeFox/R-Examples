#' Levelplot of XYZ Data
#' 
#' Plots the interpolated resistance values of the xyz 
#' data without topography (height adjustment).
#' The xyz values are obtained through inversion of the raw values.
#' @param Profile profile.
#' @param xlab label for x-axes.
#' @param ylab label for y-axes.
#' @param main title to be plotted.
#' @param col vector of colors.
#' @param breaks number of color breaks.
#' @param trafo transformation to be done on data (default: log).
#' @param backtrafo back transformation to plot correct labels (default: exp).
#' @param ... lattice levelplot arguments.
#' @export
#' @seealso \code{\link{XyzData-class}}, \code{\link{levelplotXyzHeight}},
#' \code{\link{levelplotRaw}}
#' @examples 
#' # data(sinkhole)
#' 
#' # levelplotXyz(sinkhole@profiles[[1]])
#' # levelplotLegendLabel()
#' 
#' # levelplotXyz(sinkhole@profiles[[2]])
#' # levelplotLegendLabel()
#' 
#' # levelplotXyz(sinkhole@profiles[[3]])
#' # levelplotLegendLabel()
levelplotXyz <- function(Profile, xlab="Length [m]", ylab="Depth [m]",
                         main=paste(Profile@title, "without topography"),
                         col = colors, breaks=18, trafo=log, backtrafo=exp, ...) {
  lab.breaks <- round(backtrafo(seq(trafo(min(Profile@xyzData@heightAdaption$val)),
                                    trafo(max(Profile@xyzData@heightAdaption$val)), 
                                    length.out=breaks)))
  levelplot(trafo(Profile@xyzData@seaLevel$val) ~ Profile@xyzData@seaLevel$dist * Profile@xyzData@seaLevel$depth, 
            col.regions=colorRampPalette(col),
            xlab=xlab, ylab=ylab, main=main,
            colorkey=list(
              at=as.numeric(factor(c(seq(from=0, 
                                         to=1, 
                                         length.out = breaks)))),
              labels=as.character(lab.breaks)), ...)
}

#' Levelplot of XYZ Data regarding Topography
#' 
#' Plots the interpolated resistance values of the 
#' xyz data after height adjustment.
#' The xyz values are obtained through inversion of the raw values.
#' @param xlab label for x-axes.
#' @param ylab label for y-axes.
#' @param main title to be plotted.
#' @param Profile profile.
#' @param col vector of colors.
#' @param breaks number of color breaks.
#' @param trafo transformation to be done on data (default: log).
#' @param backtrafo back transformation to plot correct labels (default: exp).
#' @param ... lattice levelplot arguments.
#' @export
#' @seealso \code{\link{XyzData-class}}, \code{\link{levelplotXyz}},
#' \code{\link{levelplotRaw}}
#' @examples 
#' # data(sinkhole)
#' 
#' # levelplotXyzHeight(sinkhole@profiles[[1]])
#' # levelplotLegendLabel()
#' 
#' # levelplotXyzHeight(sinkhole@profiles[[2]])
#' # levelplotLegendLabel()
#' 
#' # levelplotXyzHeight(sinkhole@profiles[[3]])
#' # levelplotLegendLabel()
levelplotXyzHeight <- function(Profile, xlab="Length [m]", ylab="Depth [m]",
                               main=paste(Profile@title, "with topography"),
                               col = colors, breaks=18, trafo=log, backtrafo=exp, ...) {
  lab.breaks <- round(backtrafo(seq(trafo(min(Profile@xyzData@heightAdaption$val)),
                                    trafo(max(Profile@xyzData@heightAdaption$val)), 
                                    length.out=breaks)))
  levelplot(trafo(Profile@xyzData@heightAdaption$val) ~ round(Profile@xyzData@heightAdaption$dist) * round(Profile@xyzData@heightAdaption$depth), 
            col.regions=colorRampPalette(col), 
            xlab=xlab, ylab=ylab, main=main,
            colorkey=list(
              at=as.numeric(factor(c(seq(from=0, 
                                         to=1, 
                                         length.out = breaks)))),
              labels=as.character(lab.breaks)), ...)
}

#' Levelplot Legend Label
#' 
#' Plots the label of the levelplot.
#' @param legend.lab label (default: "Resistivity").
#' @param unit unit (default: "Ohm*m").
#' @export
#' @seealso \code{\link{levelplotXyz}}, \code{\link{levelplotXyzHeight}},
#' \code{\link{levelplotRaw}}
#' @examples 
#' # data(sinkhole)
#' 
#' # levelplotRaw(sinkhole@profiles[[1]])
#' # levelplotLegendLabel()
#' 
#' # levelplotXyz(sinkhole@profiles[[2]])
#' # levelplotLegendLabel()
#' 
#' # levelplotXyzHeight(sinkhole@profiles[[3]])
#' # levelplotLegendLabel()
levelplotLegendLabel <- function(legend.lab="Resistivity", 
                                  unit=expression(paste("[", Omega, "m]"))) {
  trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
  grid.text(legend.lab, 0.2, 0, hjust=0.5, vjust=1.5)
  grid.text(unit, 0.2, 0, hjust=0.5, vjust=2.5)
  trellis.unfocus()
}