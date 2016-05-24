nullclines <- function(deriv, x.lim, y.lim, parameters = NULL, points = 101, 
                       system = "two.dim", colour = c("red", "blue"),
					   add = TRUE, ...){
  if ((is.vector(x.lim) == FALSE) | (length(x.lim) != 2)){
    stop("x.lim is not a vector of length 2 as required")
  }
  if (x.lim[2] <= x.lim[1]){
    stop("x.lim[2]  is less than or equal to x.lim[1]")
  }
  if ((is.vector(y.lim) == FALSE) | (length(y.lim) != 2)){
    stop("y.lim is not a vector of length 2 as required")
  }
  if (y.lim[2] <= y.lim[1]){
    stop("y.lim[2]  is less than or equal to y.lim[1]")
  }
  if (points <= 0){
    stop("points is less than or equal to zero")
  }
  if (!(system %in% c("one.dim", "two.dim"))){
    stop("system must either be set to one.dim or two.dim")
  }
  if (!is.vector(colour)){
    stop("colour is not a vector as required")
  }
  if (length(colour) != 2){
    if (length(colour) == 1){
      colour <- rep(colour, 2)
    }
    if (length(colour) > 2){
      colour <- colour[1:2]
    }
    print("Note: colour has been reset as required")
  }
  if (!is.logical(add)){
    stop("add must be logical")
  }
  x <- seq(from = x.lim[1], to = x.lim[2], length = points)
  y <- seq(from = y.lim[1], to = y.lim[2], length = points)
  dx <- matrix(0, ncol = points, nrow = points)
  dy <- matrix(0, ncol = points, nrow = points)
  if (system == "two.dim"){
    for (i in 1:length(x)){
      for (j in 1:length(y)){
        df <- deriv(t = 0, y = c(x[i], y[j]), parameters = parameters)
        dx[i, j] <- df[[1]][1]
        dy[i, j] <- df[[1]][2]
      }
    }
    contour(x, y, dx, levels = 0, add = add, col = colour[1], 
            drawlabels = FALSE, ...)
    contour(x, y, dy, levels = 0, add = TRUE, col = colour[2], 
            drawlabels = FALSE, ...)
  }
  if (system == "one.dim"){
    for (i in 1:length(y)){
      dy[1, i] <- deriv(t = 0, y = y[i], parameters = parameters)[[1]]
    }
    for (i in 2:length(x)){
      dy[i, ] <- dy[1, ]
    }
    contour(x, y, dy, levels = 0, add = add, col = colour[1], 
            drawlabels = FALSE, ...)
  }
  output            <- list()
  output$colour     <- colour
  output$deriv      <- deriv
  if (system == "two.dim"){
    output$dx     <- dx
  }
  output$dy         <- dy
  output$parameters <- parameters
  output$points     <- points
  output$system     <- system
  output$x.lim      <- x.lim
  output$y.lim      <- y.lim
  output$x          <- x
  output$y          <- y
  return(output)
}