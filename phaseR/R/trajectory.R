trajectory <- function(deriv, y0 = NULL, n = NULL, t.start = 0, t.end,
                       t.step = 0.01, parameters = NULL, system = "two.dim",
					   colour = "black", add = TRUE, ...){
  if (t.start < 0){
    stop("t.start is less than zero")
  }
  if (t.end <= 0){
    stop("t.end is less than or equal to zero")
  }
  if (t.start >= t.end){
    stop("t.end is less than or equal to t.start")
  }
  if (t.step <= 0){
    stop("t.step is less than or equal to zero")
  }
  if (!(system %in% c("one.dim", "two.dim"))){
    stop("system must either be set to one.dim or two.dim")
  }
  if (!is.vector(colour)){
    stop("colour is not a vector as required")
  }
  if (!is.logical(add)){
    stop("add must be logical")
  }
  if (is.null(y0) & is.null(n)){
    stop(paste("Both y0 and n cannot be NULL"))
  }
  if (!is.null(y0) & !is.null(n)){
    warning("n is non-NULL whilst y0 has also been specified")
  }
  if (is.null(y0)){
    y0 <- locator(n = n)
    if (system == "two.dim"){
      re.set <- matrix(0, ncol = 2, nrow = n)
      for (i in 1:n){
        re.set[i, ] <- c(y0$x[i], y0$y[i])
      }
      y0 <- re.set
    }
    if (system == "one.dim"){
      re.set <- rep(0, n)
      for (i in 1:n){
        re.set[i] <- y0$y[i]
      }
      y0 <- re.set
    }
  }
  if ((!is.vector(y0)) & (!is.matrix(y0))){
    stop("y0 is neither a number, vector or matrix as required")
  }
  if ((is.vector(y0)) & (length(y0) == 1) & (system == "two.dim")){
    stop("system cannot be set to two.dim, but y0 be only a number")
  }
  if ((is.matrix(y0)) & (system == "one.dim")){
    stop("system cannot be set to one.dim, and y0 be a matrix")
  }
  if (is.matrix(y0)){
    if (dim(y0)[1] > length(colour)){
      colour <- rep(colour, dim(y0)[1])
      print("Note: colour has been reset as required")
    }
    if (dim(y0)[1] < length(colour)){
      colour <- colour[1:dim(y0)[1]]
      print("Note: colour has been reset as required")
    }
  }
  if ((is.vector(y0)) & (system == "two.dim")){
    if (length(colour) != 1){
      colour <- colour[1]
      print("Note: colour has been reset as required")
    }
  }
  if ((is.vector(y0)) & (system == "one.dim")){
    if (length(y0) > length(colour)){
      colour <- rep(colour, length(y0))
      print("Note: colour has been reset as required")
    }
    if (length(y0) < length(colour)) {
      colour <- colour[1:length(y0)]
      print("Note: colour has been reset as required")
    }
  }
  t <- seq(from = t.start, to = t.end, by = t.step)
  if ((system == "two.dim") & (is.vector(y0))){
    x <- matrix(0, nrow = length(t), ncol = 1)
    y <- matrix(0, nrow = length(t), ncol = 1)
    y0 <- t(as.matrix(y0))
  }
  if ((system == "two.dim") & (is.matrix(y0))){
    x <- matrix(0, nrow = length(t), ncol = dim(y0)[1])
    y <- matrix(0, nrow = length(t), ncol = dim(y0)[1])
  }
  if (system == "one.dim"){
    y0 <- as.matrix(y0)
    x <- matrix(0, nrow = length(t), ncol = dim(y0)[1])
  }
  for (i in 1:dim(y0)[1]){
    if (system == "one.dim"){
      phase.trajectory <- ode(times = t, y = as.vector(y0[i, 1]), func = deriv,
	                          parms = parameters, method = "ode45")
    }
    if (system == "two.dim"){
      phase.trajectory <- ode(times = t, y = as.vector(y0[i, ]), func = deriv,
	                          parms = parameters, method = "ode45")
    }
    x[, i] <- phase.trajectory[, 2]
    if ((add == FALSE) & (i == 1)){
      if (system == "one.dim"){
        plot(phase.trajectory[, 1], phase.trajectory[, 2], col = colour[i],
		     type = "l", ...)
      }
      if (system == "two.dim"){
        plot(phase.trajectory[, 2], phase.trajectory[, 3], col = colour[i],
		     type = "l", ...)
        y[, i] <- phase.trajectory[, 3]
      }
    }
    if ((add == TRUE) | (i > 1)){
      if (system == "one.dim"){
        lines(phase.trajectory[, 1], phase.trajectory[, 2], col = colour[i],
		      type = "l", ...)
      }
      if (system == "two.dim"){
        lines(phase.trajectory[, 2], phase.trajectory[, 3], col = colour[i],
		      type = "l", ...)
        y[, i] <- phase.trajectory[, 3]
      }
    }
  }
  if (system == "two.dim"){
    points(y0[, 1], y0[, 2], col = colour, ...)
  }
  if (system == "one.dim"){
    points(rep(t.start, dim(y0)[1]), y0[, 1], col = colour, ...)
  }
  output            <- list()
  output$colour     <- colour
  output$deriv      <- deriv
  output$parameters <- parameters
  output$system     <- system
  output$t.start    <- t.start
  output$t.step     <- t.step
  output$t.end      <- t.end
  output$t          <- t
  if (system == "two.dim"){
    output$x        <- x
  }
  output$y0         <- y0
  if (system == "one.dim"){
    output$y        <- x
  }
  if (system == "two.dim"){
    output$y        <- y
  }
  return(output)
}