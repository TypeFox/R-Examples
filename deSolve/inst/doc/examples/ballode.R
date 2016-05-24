## =============================================================================
## A bouncing ball; ode with event location
## =============================================================================
require(deSolve)
#-----------------------------
# the model function
#-----------------------------
ballode<- function(t, y, parms) {
  dy1 <- y[2]
  dy2 <- -9.8
  list(c(dy1, dy2))
}

#-----------------------------
# the root and event function
#-----------------------------
# event triggered when the ball hits the ground (height = 0)
root <- function(t, y, parms) y[1]

# bouncing
event <- function(t, y, parms) {
  y[1] <- 0
  y[2] <- -0.9 * y[2]
 return(y)
}

#-----------------------------
# initial values and times
#-----------------------------
yini  <- c(height = 0, v = 20)
times <- seq(0, 40, 0.01)

#-----------------------------
# solve the model
#-----------------------------
out   <- lsodar(times = times, y = yini, func = ballode, parms = NULL,
  events = list(func = event, root = TRUE), rootfun = root)

out2   <- radau(times = times, y = yini, func = ballode, parms = NULL,
  events = list(func = event, root = TRUE), rootfun = root)  #         , verbose=TRUE

attributes(out)$troot
attributes(out2)$troot
#-----------------------------
# display, plot results
#-----------------------------
plot(out, which = "height", type = "l", lwd = 2, main = "bouncing ball", ylab = "height")

