

#' Plot conditional effects
#' 
#' Can be used to make a conditional effects plot with an effect function on the
#' y axis and a covariate on the x axis. \code{ggplot2} is used to create the plot.
#' 
#' @param obj Object of class \code{effectlite} obtained from fitting an effect 
#' model using \code{\link[EffectLiteR]{effectLite}} 
#' @param zsel Name of a covariate (character string) plotted on the x-axis.
#' @param gxsel Name of an effect function (character string) plotted on the y-axis.
#' @param colour Name of a covariate (character string) used as colour variable 
#' in the plot.
#' @return Object of class \code{c("gg", "ggplot")}.
#' @examples
#' m1 <- effectLite(y="dv", x="x", k="k1", z="z1", control="control", data=example01)
#' conditionalEffectsPlot(m1, zsel="z1", gxsel="g1", colour="k1")
#' 
#' @export
conditionalEffectsPlot <- function(obj, zsel, gxsel="g1", colour=""){
  
  stopifnot(class(obj) == "effectlite")  
  
  condeffects <- obj@results@condeffects
  
  stopifnot(zsel %in% names(condeffects))
  stopifnot(gxsel %in% names(condeffects))
  
  yselected <- round(condeffects[[gxsel]],4)    
  zselected <- condeffects[[zsel]]
  colourselected <- condeffects[[colour]]
  
  g1label <- "(K,Z)"
  if(!"K" %in% names(condeffects)){g1label <- "(Z)"}
  
  p <- ggplot2::qplot(y=yselected, x=zselected, 
                      data=condeffects,
                      ylab=paste0(gxsel,g1label),
                      xlab=zsel,                 
                      main=paste0("Estimated regression of ",
                                  paste0(gxsel,g1label), " on ", 
                                  zsel))
  p <- p + ggplot2::geom_smooth(method="loess")
  if(!is.null(colourselected)){
    p <- p + ggplot2::geom_point(ggplot2::aes(colour=colourselected))
    p <- p + ggplot2::guides(colour = ggplot2::guide_legend(colour))
  }
  p <- p + ggplot2::theme_bw()
  
  return(p)
  
}

