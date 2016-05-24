########################################################################
## Random walk model with barriers
## (C) Thomas Petzoldt, License GPL >= 2
########################################################################

library(simecol)

dm <- rwalkModel(
  main = function(time, init, parms, inputs) {
    speed   <- parms$speed
    xleft   <- parms$area[1]
    xright  <- parms$area[2]
    ybottom <- parms$area[3]
    ytop    <- parms$area[4]

    ## change barriers at time > tcrit
    if (time > parms$tcrit)
      barriers <- inputs$barriers1
    else
      barriers <- inputs$barriers2

    x <- init$x  # x coordinate
    y <- init$y  # y coordinate
    n <- length(x)

    ## Random Walk (rectangular, one cell per time step only)
    dx <- sample(c(-1, 0, 1), n, replace = TRUE)
    dy <- sample(c(-1, 0, 1), n, replace = TRUE)
    xnew  <- x + dx
    ynew  <- y + dy

    ## Find position of individuals in "environment matrix"  (here: barriers)
    i.j <- array(c(pmax(1, xnew), pmax(1, ynew)), dim = c(n, 2))

    ## if new position is a barrier, then stay at old location,
    ## otherwise assign new location
    ## (it would be also possible to handle this for x and y separately)
    x <- ifelse(barriers[i.j] == 1, x, xnew)
    y <- ifelse(barriers[i.j] == 1, y, ynew)

    ## return new state (data frame containing the individuals)
    data.frame(x = x, y = y, ID = init$ID, genotype = init$genotype)
  },
  times  = c(from = 0, to = 300, by = 1),
  parms  = list(ninds = 50, area = c(1, 100, 1, 100), tcrit=200),
  solver = "iteration",
  initfunc = function(obj) {
    ## initfunc creates the start population and the environment matrices (barriers)
    p <- parms(obj)

    population <- with(p,
      data.frame(x = sample((area[1]+1):(area[2]-1), ninds, replace = TRUE),
                 y = sample((area[3]+1):(area[4]-1), ninds, replace = TRUE),
                 ID = 1:ninds)
    )

    # lower population has genotype "1", upper population has genotype "2"
    population$genotype <- ifelse(population$y < 50, 1, 2)

    init(obj) <- population

    barriers1 <- matrix(0, nrow = 100, ncol = 100)
    barriers1[p$area[1], ] <- 1  # left boundary
    barriers1[p$area[2], ] <- 1  # right boundary
    barriers1[, p$area[3]] <- 1  # bottom boundary
    barriers1[, p$area[4]] <- 1  # top boundary

    barriers2 <- barriers1       # make a copy
    barriers2[c(1:40, 60:100),  50] <- 1  # "barriers2" segments populations

    ## assign two barrier matrices, one for the first phase, one for the second
    inputs(obj) <- list(barriers1 = barriers1, barriers2 = barriers2)
    obj
  }
)

observer(dm) <- function(state, time, i, out, obj) {
  ## numerical output to the screen
  cat("time =", time, "\n")
  ## animation
  p   <- parms(obj)
  inp <- inputs(obj)
  with(p, {
    ## switch between barriers at time > tcrit
    if (time > tcrit)
      barriers <- inp$barriers1
    else
      barriers <- inp$barriers2

    image(barriers, col=c("wheat", "grey"), axes = FALSE)
    points((state$x - area[1])/(area[2] - area[1]),
           (state$y - area[3])/(area[4] - area[3]),
    xlab="x", ylab="y", pch=16, col=c("red", "blue")[state$genotype])
    box() # cosmetics
  })
  ## return the state --> iteration stores it in "out"
  state
}


parms(dm)["ninds"] <- 200

par(mar=rep(1, 4))
dm <- sim(dm)


