ModelPlot = function (model=NULL, user.func = NULL, dimensions = list(x1 = NULL, x2 = NULL, 
                     x3 = NULL), slice = NULL, lims = rep(0, 6), constraints = FALSE, 
                      constraint.pars = list(lty = 2, lwd = 2), contour = FALSE, 
                      contour.pars = list(lwd = 0.5, cex.lab = 1.3), cuts = 10, 
                      at = NULL, res = 300, pseudo = FALSE, fill = FALSE, color.palette = heat.colors, 
                      main = NULL, axislabs = c("Fraction X1", "Fraction X2", "Fraction X3"), 
                      axislab.pars = list(), axislab.offset = 0, cornerlabs = c("X1", 
                      "X2", "X3"), cornerlab.pars = list(), grid = TRUE, grid.pars = list(col = "darkgrey", 
                      lty = 3, lwd = 0.5), colorkey = FALSE, labels = TRUE, 
                      label.style = "align", ...) 
{
  if (!is.null(slice) & (is.null(slice$mix.vars) & is.null(slice$process.vars))) {
    stop("Slice argument must be of form: slice = list(mix.vars = c(...), process.vars = c(...))")
  }
  if (is.null(slice$mix.vars)) {
    mx = 1
  } else {
    mx = 1 - sum(unlist(slice$mix.vars))
  }
  
  if ((constraints | pseudo) & sum(lims == rep(0, 6)) == 6) {
    stop("Component limits must be specified with the 'lims' argument to plot constraints or to use pseudo components")
  }
  if ((constraints & (lims[1] > mx | lims[3] > mx | lims[5] > 
                        mx))) {
    stop("Lower constraints must be less than the sum of fixed component proportions.  Constraints are inconsistent")
  }
  if(!is.null(model) && is.null(user.func)){
    user.func = function(grid){
      ygrid = predict(model,grid)
      return(ygrid)
    }
  }
  if (is.null(user.func)) 
    stop("ModelPlot requires either the model argument or the user.func argument")
  if (length(dimensions) != 3) 
    stop("Ternary plot dimensions must be identified")
  if (is.logical(colorkey)) {
    if (colorkey) {
      colorkey = list(axis.line = list(col = 1), axis.text = list(col = 1))
    }
  }
  else {
    if (is.list(colorkey)) {
      if (is.null(colorkey$axis.line)) {
        colorkey$axis.line = list(col = 1)
      }
      if (is.null(colorkey$axis.text)) {
        colorkey$axis.text = list(col = 1)
      }
    }
    else {
      warning("colorkey argument must be a logical or a list.  Argument will be ignored.")
      colorkey = FALSE
    }
  }
  base = high = NULL
  l.bnds <- lims[seq(1, 5, by = 2)]
  if (sum(l.bnds >= mx) > 0) 
    stop("The lower bound on at least one component on the ternary plot is larger than the proportion of the mixture not occupied by the fixed component(s).  The constraints are not consistent.")
  trian <- expand.grid(base = seq(0, 1, l = res), high = seq(0, 
                                                             sin(pi/3), l = res))
  trian <- subset(trian, (base * sin(pi/3) * 2) >= high)
  trian <- subset(trian, ((1 - base) * sin(pi/3) * 2) >= high)
  new2 <- data.frame(x1 = trian$high * 2/sqrt(3))
  new2$x3 <- trian$base - trian$high/sqrt(3)
  new2$x2 <- 1 - new2$x1 - new2$x3
  new2 = new2[, c(1, 3, 2)]
  if (!pseudo) {
    l.bnds2 = rep(0, 3)
  }
  else {
    l.bnds2 = l.bnds
  }
  sum.bnds <- sum(l.bnds2)
  new2$x3 <- l.bnds2[3] + (mx - sum.bnds) * new2$x3
  new2$x2 <- l.bnds2[2] + (mx - sum.bnds) * new2$x2
  new2$x1 <- l.bnds2[1] + (mx - sum.bnds) * new2$x1
  names(new2) = unlist(dimensions)
  
  if (!is.null(slice$mix.vars)) {
    for (i in 1:length(slice$mix.vars)) {
      eval(parse(text = paste("new2 = cbind(new2,",names(slice$mix.vars[i]),"=slice$mix.vars[i],row.names=NULL)")))
    }
  }
  if (!is.null(slice$process.vars)) {
    for (i in 1:length(slice$process.vars)) {
      eval(parse(text = paste("new2 = cbind(new2,",names(slice$process.vars[i]),"=slice$process.vars[i],row.names=NULL)")))
    }
  }
  
  trian$w = user.func(grid = new2, ...)
  for (i in seq(2, 6, by = 2)) {
    if (lims[i] == 0) 
      lims[i] = 1
  }
  grade.trellis <- function(...) {
    from = 0.2
    to = 0.8
    step = 0.2
    x1 <- seq(from, to, step)
    x2 <- x1/2
    y2 <- x1 * sqrt(3)/2
    x3 <- (1 - x1) * 0.5 + x1
    y3 <- sqrt(3)/2 - x1 * sqrt(3)/2
    labx1 <- l.bnds2[3] + (mx - sum.bnds) * x1
    labx2 <- l.bnds2[2] + (mx - sum.bnds) * x1
    labx3 <- l.bnds2[1] + (mx - sum.bnds) * x1
    panel.segments(x1, 0, x2, y2, ...)
    panel.text(x1, 0, label = labx1, pos = 1)
    panel.segments(x1, 0, x3, y3, ...)
    panel.text(x2, y2, label = rev(labx2), pos = 2)
    panel.segments(x2, y2, 1 - x2, y2, ...)
    panel.text(x3, y3, label = rev(labx3), pos = 4)
  }
  plot.constraints <- function(...) {
    if (constraints) {
      x1 <- (lims[2] - l.bnds2[1])/(mx - sum.bnds)
      x2 <- x1/2
      y2 <- x1 * sqrt(3)/2
      if (x1 < 1) 
        panel.segments(x2, y2, 1 - x2, y2, ...)
      x1 <- 1 - (lims[4] - l.bnds2[2])/(mx - sum.bnds)
      x2 <- x1/2
      y2 <- x1 * sqrt(3)/2
      if (x1 > 0) 
        panel.segments(x1, 0, x2, y2, ...)
      x1 <- (lims[6] - l.bnds2[3])/(mx - sum.bnds)
      x2 <- (1 - x1)/2 + x1
      y2 <- sqrt(3)/2 - x1 * sqrt(3)/2
      if (x1 < 1) 
        panel.segments(x1, 0, x2, y2, ...)
      x1 <- (lims[1] - l.bnds2[1])/(mx - sum.bnds)
      x2 <- x1/2
      y2 <- x1 * sqrt(3)/2
      if (x1 > 0) 
        panel.segments(x2, y2, 1 - x2, y2, ...)
      x1 <- 1 - (lims[3] - l.bnds2[2])/(mx - sum.bnds)
      x2 <- x1/2
      y2 <- x1 * sqrt(3)/2
      if (x1 < 1) 
        panel.segments(x1, 0, x2, y2, ...)
      x1 <- (lims[5] - l.bnds2[3])/(mx - sum.bnds)
      x2 <- (1 - x1) * 0.5 + x1
      y2 <- sqrt(3)/2 - x1 * sqrt(3)/2
      if (x1 > 0) 
        panel.segments(x1, 0, x2, y2, ...)
    }
  }
  if (is.null(at)) {
    levelplot.main = list(x = w ~ base * high, data = trian, 
                          aspect = "iso", xlim = c(-0.1, 1.1), ylim = c(-0.1, 
                                                                        0.96), main = main, xlab = NULL, ylab = NULL, 
                          contour = contour, cuts = cuts, labels = labels, 
                          label.style = label.style, pretty = TRUE, region = fill, 
                          col.regions = color.palette(n = cuts + 1), par.settings = list(axis.line = list(col = NA), 
                                                                                         axis.text = list(col = NA)), colorkey = colorkey)
  }
  else {
    levelplot.main = list(x = w ~ base * high, data = trian, 
                          aspect = "iso", xlim = c(-0.1, 1.1), ylim = c(-0.1, 
                                                                        0.96), main = main, xlab = NULL, ylab = NULL, 
                          contour = contour, at = at, labels = labels, label.style = label.style, 
                          pretty = TRUE, region = fill, col.regions = color.palette(n = length(at) + 
                                                                                      1), par.settings = list(axis.line = list(col = NA), 
                                                                                                              axis.text = list(col = NA)), colorkey = colorkey)
  }
  levelplot.args1 = c(levelplot.main, contour.pars, panel = function(...) {
    panel.levelplot(...)
    panel.segments(c(0, 0, 0.5), c(0, 0, sqrt(3)/2), c(1, 
                                                       1/2, 1), c(0, sqrt(3)/2, 0), lwd = 2)
    if (grid) {
      do.call(grade.trellis, grid.pars)
    }
    do.call(plot.constraints, constraint.pars)
    do.call(panel.text, c(0, 0, label = cornerlabs[2], pos = 2, 
                          cornerlab.pars))
    do.call(panel.text, c(1/2, sqrt(3)/2, label = cornerlabs[1], 
                          pos = 3, cornerlab.pars))
    do.call(panel.text, c(1, 0, label = cornerlabs[3], pos = 4, 
                          cornerlab.pars))
    do.call(panel.text, c(0.5, -0.075, axislabs[3], axislab.pars))
    do.call(panel.text, c(0.15 - axislab.offset, 0.5, axislabs[2], 
                          srt = 60, axislab.pars))
    do.call(panel.text, c(0.85 + axislab.offset, 0.5, axislabs[1], 
                          srt = -60, axislab.pars))
  })
  p = do.call(levelplot, levelplot.args1)
  print(p)
}