numericalSolution <- function(deriv, y0 = NULL, t.start = 0, t.end,
                              t.step = 0.01, parameters = NULL, type = "two",
							  colour = rep("black", 2), grid = TRUE, ...){
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
  if (is.null(y0)){
    y0 <- locator(n = 1)
  }
  if (!is.vector(y0)){
    stop("y0 is not a vector as required")
  }
  if (length(y0) == 1){
    stop("y0 should be a vector of length two")
  }
  if (!(type %in% c("one", "two"))){
    stop("type must either be set to one or two")
  }
  if (!is.vector(colour)){
    stop("colour is not a vector as required")
  }
  if (length(colour) == 1){
    colour <- rep(colour, length(y0))
    print("Note: colour has been reset as required")
  }
  if (length(colour) > 2){
    colour <- colour[1:2]
    print("Note: colour has been reset as required")
  }
  if (!is.logical(grid)){
    stop(paste("grid must either be set to TRUE or FALSE"))
  }
  t <- seq(from = t.start, to = t.end, by = t.step)
  phase.trajectory <- ode(times = t, y = y0, func = deriv, 
                          parms = parameters, method = "ode45")
  x <- phase.trajectory[, 2]
  y <- phase.trajectory[, 3]
  if (type == "one"){
    plot(t, x, col = colour[1], type = "l", ...)
    lines(t, y, col = colour[2], type = "l", ...)
    if (grid == TRUE){
      grid()
    }
  }
  if (type == "two"){
    old.par <- par(no.readonly = TRUE)
    par(mfcol = c(2, 1), oma = c(0, 0, 2, 0))
    plot(t, x, col = colour[1], type = "l", ...)
    if (grid == TRUE){
      grid()
    }
    plot(t, y, col = colour[2], type = "l", ...)
    if (grid == TRUE){
      grid()
    }
    par <- old.par
  }
  output            <- list()
  output$colour     <- colour
  output$parameters <- parameters
  output$t.start    <- t.start
  output$t.step     <- t.step
  output$t.end      <- t.end
  output$t          <- t
  output$type       <- type
  output$x          <- x
  output$y0         <- y0
  output$y          <- y
  return(output)
}