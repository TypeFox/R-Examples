eplot <-
function(data, element, panel=NULL, groups=NULL, ...) {
  mybg <- trellis.par.get('strip.background')
  mybg$col <- 'white'
  trellis.par.set('strip.background', mybg)
  class(data) <- 'data.frame'
  tmp.depth <- depths(data$depth, data$name)[,1]
  variable <- data[,which(names(data)==element)]
  if (is.null(panel)) {panel.variable <- data$Profile} else {panel.variable <- panel}
  if (is.null(groups)==TRUE) {
    xyplot(tmp.depth~variable|panel.variable, ...)
  } else {
    xyplot(tmp.depth~variable|panel.variable, groups=groups, ...)
  }
}
