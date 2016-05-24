phasePortrait <- function(deriv, y.lim, y.step = 0.01, parameters = NULL,
                          points = 10, frac = 0.75, arrow.head = 0.075,
						  colour = "black", xlab = "y", ylab = "f(y)", ...){
    if ((is.vector(y.lim) == FALSE) | (length(y.lim) != 2)){
        stop("y.lim is not a vector of length 2 as required")
    }
    if (y.lim[2] <= y.lim[1]){
        stop("y.lim[2]  is less than or equal to y.lim[1]")
    }
    if (y.step <= 0){
        stop("y.step is less than or equal to zero")
    }
    if (!is.vector(colour)){
        stop("colour is not a vector as required")
    }
    if (length(colour) > 1){
        colour <- colour[1]
        print("Note: colour has been reset as required")
    }
    y <- seq(from = y.lim[1], to = y.lim[2], by = y.step)
    dy <- rep(0, length(y))
    for (i in 1:length(y)){
        dy[i] <- deriv(t = 0, y = y[i], parameters = parameters)[[1]]
    }
    plot(y, dy, col = colour, type = "l", xlab = xlab, ylab = ylab, ...)
    y.arrows <- seq(from = y.lim[1], to = y.lim[2], length = points)
    dy.arrows <- rep(0, points)
    y.shift <- 0.5*frac*(y.arrows[2] - y.arrows[1])
    for (i in 1:points){
        dy.arrows[i] <- deriv(t = 0, y = y.arrows[i],
		                      parameters = parameters)[[1]]
    }
    pos <- which(dy.arrows > 0)
    arrows(y.arrows[pos] - y.shift, rep(0, length(y.arrows[pos])), 
           y.arrows[pos] + y.shift, rep(0, length(y.arrows[pos])), 
           length = arrow.head, col = colour, ...)
    neg <- which(dy.arrows < 0)
    arrows(y.arrows[neg] + y.shift, rep(0, length(y.arrows[neg])), 
           y.arrows[neg] - y.shift, rep(0, length(y.arrows[neg])), 
           length = arrow.head, col = colour, ...)
    output            <- list()
    output$arrow.head <- arrow.head
    output$colour     <- colour
    output$deriv      <- deriv
    output$dy         <- dy
    output$frac       <- frac
    output$parameters <- parameters
    output$y.step     <- y.step
    output$y.lim      <- y.lim
    output$y          <- y
	return(output)
}