plot.clusterTree <-
function(x, type = "lambda", color = NULL, add = FALSE, ...) {
  
  if (class(x) != "clusterTree") {
    stop("x should be an object of class clusterTree")
  }
  if (type != "lambda" && type != "r" && type != "kappa" && type != "alpha") {
    stop("type should be either 'r', 'lambda', 'alpha' or 'kappa'")
  }
  if (!is.logical(add)) {
    stop("add should be logical")
  }
  
  id <- x[["id"]]
  base <- x[["Xbase"]]
  
  if (type == "lambda"){
    bottom <- x[["lambdaBottom"]]
    top <- x[["lambdaTop"]]
  } else if (type == "r") {
    bottom <- -x[["rBottom"]]
    top <- -x[["rTop"]]
  } else if (type == "kappa") {
    bottom <- x[["kappaBottom"]]
    top <- x[["kappaTop"]]
  } else if (type == "alpha") {
    bottom <- -x[["alphaBottom"]]
    top <- -x[["alphaTop"]]
  }

  Range <- range(c(bottom, top))
  

  AT <- c(bottom, top)
  labels <- -round(AT, 3)
  if (type == "lambda" || type == "kappa") {
    labels <- round(AT,3)
  }
  
  sons <- x[["sons"]]
  
  Ylim <- c(min(bottom[which(!is.na(bottom))]), max(top[which(!is.na(top))]))
  
  if (!add) {
    graphics::plot(c(0, 0), c(Ylim[1], Ylim[2]), type = "n", xlim = c(0, 1),
        ylab = type, xlab = "", axes = FALSE, ...)
    graphics::axis(1)
    graphics::axis(2, at = AT, labels = labels)
  }
  #vertical lines
  if (is.null(color)) {
    color <- id
  }
  graphics::segments(base, bottom, base, top, col = color, lwd = 3)

  ## now the horizontal lines
  if (length(sons) > 0) { 
    for (i in seq(along = sons)) {
      if (!is.null(sons[[i]]) && length(sons[[i]]>0) && !is.na(sons[[i]])) {
        x <- c(min(base[sons[[i]]]) , max(base[sons[[i]]]))
        y <- bottom[sons[[i]][1]]
        graphics::segments(x[1], y, x[2], y, lwd = 3,
            col = ifelse(is.null(color), 1, color))
      } 
    }
  }
}
