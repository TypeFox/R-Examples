########################################################################
## Random walk model, see also diffusion_B which uses "initfunc"
########################################################################

diffusion <- new("rwalkModel",
  main = function(time, init, parms, inputs) {
    # inputs  <- obj@inputs
    speed   <- parms$speed
    xleft   <- parms$area[1]
    xright  <- parms$area[2]
    ybottom <- parms$area[3]
    ytop    <- parms$area[4]

    x <- init$x  # x coordinate
    y <- init$y  # y coordinate
    a <- init$a  # angle (in radians)
    n <- length(a)

    ## Rule 1: respect environment (grid as given in "inputs")
    ## 1a) identify location on "environmental 2D grid" for each individual
    i.j <- array(c(pmax(1, ceiling(x)), pmax(1, ceiling(y))), dim=c(n, 2))
    # print(i.j)
    ## 1b) speed dependend on "environmental conditions"
    speed <- speed * inputs[i.j]
    ## Rule 2: Random Walk
    a  <- (a + 2 * pi / runif(a)) %% (2 * pi)
    dx <- speed * cos(a)
    dy <- speed * sin(a)
    x  <- x + dx
    y  <- y + dy
    ## Rule 3: Wrap Around
    x <- ifelse(x > xright, xleft, x)
    y <- ifelse(y > ytop, ybottom, y)
    x <- ifelse(x < xleft, xright, x)
    y <- ifelse(y < ybottom, ytop, y)
    data.frame(x=x, y=y, a=a)
  },
  times  = c(from=0, to=100, by=1),
  parms  = list(ninds=50, speed = 1, area = c(0, 100, 0, 100)),
  solver = "iteration"
)

diffusionInit <- function(obj) {
  ninds   <- obj@parms$ninds
  xleft   <- obj@parms$area[1]
  xright  <- obj@parms$area[2]
  ybottom <- obj@parms$area[3]
  ytop    <- obj@parms$area[4]
  data.frame(x = runif(ninds) * (xright - xleft) + xleft,
             y = runif(ninds) * (ytop - ybottom) + ybottom,
             a = runif(ninds) * 2 * pi)
}

diffusionInputs <- function() {
  inp <- matrix(1, nrow=100, ncol=100)
  inp[, 45:55] <- 0.2
  inp
}

diffusion@init   <- diffusionInit(diffusion)
diffusion@inputs <- diffusionInputs()

