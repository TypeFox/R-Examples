########################################################################
## Stochastic Cellular Automaton with Weighted Neighborhood Matrix
########################################################################

library(simecol)
CA <- gridModel(
    main = function(time, init, parms) {
      z <- init
      n <- nrow(z)
      m <- ncol(z)
      ret <- with(parms,{
        ## rule 1: reproduction
        ## 1.1 which cells are adult? (only adults can generate)
        ad <- ifelse(z >= adult & z < old, z, 0)
        ## 1.2 how much (weighted) adult neighbours has each cell?
        nb <- neighbours(ad, wdist = wdist)
        ## 1.3 a proportion of the seeds develops juveniles
        ## simplified version, you can also use probabilities
        genprob <- nb * runif(nb) * ci
        zgen  <- ifelse(z==0 & genprob >= 1, 1, 0)
        ## rule 2: growth and survival of juveniles
        zsurvj <- ifelse(z >= 1     & z < adult & runif(z) <= pj, z+1, 0)
        ## rule 2: growth and survival of adults
        zsurva <- ifelse(z >= adult & z < old   & runif(z) <= pa, z+1, 0)
        ## rule 2: growth and survival of senescent
        zsurvs <- ifelse(z >= old               & runif(z) <= ps, z+1, 0)
      
        ## make resulting grid of complete population
        z     <- zgen + zsurvj + zsurva + zsurvs
        if (max(z)==0) stop("extinction", call.=FALSE)
        return(z)
      })
    },
    parms = list(
      wdist = matrix(c(0.5,0.5,0.5,0.5,0.5,
                  0.5,1.0,1.0,1.0,0.5,
                  0.5,1.0,1.0,1.0,0.5,
                  0.5,1.0,1.0,1.0,0.5,
                  0.5,0.5,0.5,0.5,0.5),nrow=5),
      pj = 0.99,  # survival probability of juveniles
      pa = 0.99,  # survival probability of adults
      ps = 0.1,   # survival probability of senescent
      ci = 1.0,   # "seeding constant" 
      adult = 5,  # age of adolescence
      old   = 10  # age of senescence
    ),
    times = c(from=1, to=10, by=1),
    init = matrix(0, nrow=80, ncol=80),
    solver = "iteration",
    initfunc = function(obj) {
      init(obj)[38:42,38:42] <- 5 # deterministic seed in the middle of the grid
      obj
  }
)


mycolors <- function(n) {
  col <- c("wheat", "darkgreen")
  if (n>2) col <- c(col, heat.colors(n-2))
  col
}

CA.save <- CA

times(CA) <- c(to=100)
CA <- sim(CA)         # takes some time

lastZ <- out(CA)[[50]]
maxz  <- max(lastZ)

plot(CA, delay=50, index=90, col=mycolors(maxz+1), axes=F)

# image(out(CA)[[2]])
# parms(CA) <- c(ps = 0.99)
# parms(CA) <- c(pj = 0.99)
# parms(CA) <- c(old= 20)

## alternative weight matrix for neighbourhood determination
## (square matrix with uneven number of rows and columns)
maxdist <- 20
wdist <- matrix(0, nrow = 2 * maxdist + 1,
                   ncol = 2 * maxdist + 1)
for (i in -maxdist:maxdist) {
  for (j in -maxdist:maxdist) {
    # inverse euclidean distance from center
    wdist[i + maxdist + 1, j + maxdist + 1] <- 1/sqrt(i^2 + j^2)
    # exponential decrease
    #wdist[i + maxdist + 1, j + maxdist + 1] <- exp(-sqrt(i^2 + j^2))
  }
}

wdist[maxdist + 1, maxdist + 1] <- 0 # middle cell
image(wdist)

parms(CA)$wdist <- wdist

CA <- sim(CA)
o <- out(CA)

abundance <- sapply(o, function(x) sum(x > 0))
age       <- sapply(o, function(x) sum(x)/sum(x > 0))

opar <- par(mfrow = c(2, 2))
plot(abundance)
plot(abundance, log = "y")
plot(age, ylab = "mean age")
hist(o[[50]], main="age distribution")

par(opar)
plot(CA, delay = 50, index = 50, col = mycolors(maxz + 1))
