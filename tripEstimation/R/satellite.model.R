"satellite.model" <-
function(day, X,
         proposal.x,proposal.z,
         mask.x,mask.z,
         fix.release=TRUE,fix.recapture=TRUE,
         start.x,start.z,
         posn.sigma=1,
         behav.dist="gamma", behav.mean, behav.sd,
         proj.string = "+proj=longlat"
         ) {

if (length(day) != nrow(X)) stop("length of times should match number of locations")

  ##
  ## Utility functions
  ##

  ## (Great Circle) Distance
  dist <- function(a,b) {
    r <- cos(pi/180*a[,2])*cos(pi/180*b[,2])*
      cos(pi/180*(b[,1]-a[,1]))+sin(pi/180*a[,2])*sin(pi/180*b[,2])
    6378.137*acos(pmin.int(r, 1))
  }

  # no transformation
  transf <- function(x, inv = FALSE) x


  if (!any(grep("longlat", proj.string))) {
      ## removed by namespace MDS 2011-10-06
      ## require(rgdal)

    if (!any(grep("km", proj.string))) {
      warning("Distances will be calculated in the units of the coordinate system, but kilometres are assumed for the units of speed calculation.")
      print(CRSargs(CRS(proj.string)))
    }
    dist <- function(a, b) {
      sqrt(rowSums((a - b)^2)) }

    transf <-  function(x, inv = FALSE) {
      rubbish <- capture.output(res <- project(x, proj.string, inv))
      res
    }
  }
  ##
  ## Positional contribution to the log posterior
  ##
  logp.position <- function(x) {
    x[,1:2] <- transf(x[,1:2], inv = TRUE)
    rowSums(dnorm(X, x, posn.sigma, log=TRUE))
  }


  ##
  ## Behavioural contribution to the log posterior
  ##

  ## Times between observations
  dt <- diff(unclass(day)/3600)

  if(behav.dist=="gamma") {
    alpha <- behav.mean^2/behav.sd^2
    beta <- behav.mean/behav.sd^2
    logp.behavioural <- function(k,x1,z,x2) {
      ## Average speed from x1 to z to x2
      spd <- pmax.int(dist(x1,z)+dist(z,x2), 1e-06)/dt[k]
      dgamma(spd, alpha, beta, log = TRUE)
    }
  } else {
    log.sigma <- sqrt(log(1+behav.mean^2/behav.sd^2))
    log.mu <- log(behav.mean)-log.sigma^2/2
    logp.behavioural <- function(k, x1, z, x2) {
      ## Average speed from x1 to z to x2
      spd <- (dist(x1, z) + dist(z, x2))/dt[k]
      dnorm(spd, log.mu, log.sigma, log = TRUE)
    }
  }

  ##
  ## Locations to be held fixed
  ##
  n <- length(day)
  fixed.x <- rep(FALSE,n)
  fixed.x[1] <- fix.release
  fixed.x[n] <- fix.recapture


  ## Return a list with all model components
  list(## Number of locations
       n = n,
       ## The function for generating proposal x's
       proposal.x=proposal.x,
       ## The function for generating proposal z's
       proposal.z=proposal.z,
       ## The mask for the x's
       mask.x=mask.x,
       ## The mask for the z's
       mask.z=mask.z,
       ## Positional contribution to the log posterior
       logp.position=logp.position,
       ## Behavioural contribution to the log posterior
       logp.behavioural=logp.behavioural,
       ## Locations to be held fixed
       fixed.x=fixed.x,
       ## Suggested starting points
       start.x=start.x,
       start.z=start.z)
}

