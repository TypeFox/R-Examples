"solar.model" <-
function(segments,day,light,
	proposal.x, proposal.z,
         mask.x,mask.z,
         fix.release=TRUE,fix.recapture=TRUE,
         calibration,
         light.sigma=7,k.sigma=10,
         behav = "speed",
         behav.dist="gamma", behav.mean, behav.sd,
         proj.string = "+proj=longlat",
         ekstrom = c(-5, 3, light.sigma),
         ekstrom.limit = "light"
         ) {

## add fakes here
if (fix.release) {
	segments <- c(segments[1] - 1, segments)
	day <- c(day[1] - 12*3600, day)
	light <- c(0, light)
}

if (fix.recapture) {
	segments <- c(segments, segments[length(segments)] + 1)
	day <- c(day, day[length(day)] + 12*3600)
	light <- c(light, 0)
}





## check the data
nn <- c(length(segments), length(day), length(light))
if (any(diff(nn) != 0)) stop("data of unequal lengths")
if (nn[1] < 100) stop(paste("only", nn[1], "light levels!  Check the input data. "))
test.c <- calibration(-10:10)
if (length(test.c) < 10) stop("less than 10 values from calibration, probably not enough")
if (any(c(light.sigma, k.sigma, behav.mean, behav.sd) <= 0 )) stop("variance/mean of zero or less")

## do a check on the proposals

dmx <- dim(proposal.x(matrix(1, 10)))
dmz <- dim(proposal.z(matrix(1, 10)))

if (dmx[1] != length(unique(segments))) stop("number of segments do not match X-proposals")
if (dmz[1] != (dmx[1] - 1)) stop("number of X-proposals should be one more than Z-proposals")

  ## assume long/lat

  ## (Great Circle) Distance in kilometres
  dist <- function(a,b) {
    r <- cos(pi/180*a[,2])*cos(pi/180*b[,2])*
      cos(pi/180*(b[,1]-a[,1]))+sin(pi/180*a[,2])*sin(pi/180*b[,2])
    6378.137*acos(pmin.int(r, 1))
  }

  # no transformation
  transf <- function(x, inv = FALSE) x


  if (!any(grep("longlat", proj.string))) {
##    require(rgdal) removed by namespace addition 2011-10-09 MDS

    if (!any(grep("km", proj.string))) {
      warning("Distances will be calculated in the units of the coordinate system, but kilometres are assumed for speed.")
      print(CRSargs(CRS(proj.string)))
    }
    dist <- function(a, b) {
      sqrt(rowSums((a - b)^2)) }

    transf <-  function(x, inv = FALSE) {
      rubbish <- capture.output(res <- project(x, proj.string, inv))
      res
    }
  }



  ## we work in hours
  dt.scale <- 3600

  ## Predicted light levels
  light.predict <- function(x,day) {
    ## Calculate elevations - note we must expand locations per
    ## segment to locations per observation
   elev <- elevation(x[,1],x[,2],solar(day))

    ## Expected (log scale) light level given tag calibration and the
    ## attenuation offset
   calibration(elev)+x[,3]
 }


  ##
  ## Positional contribution to the log posterior
  ##

  ## The vector of indices is maps segments to observations.
  segments <- factor(segments)
  is <- unclass(segments)

  ## Precompute solar constants
  sun <- solar(day)

  light.sigma <- rep(light.sigma, length(segments))
  if (ekstrom[3] != light.sigma[1]) {

  	if (ekstrom.limit == "elevation") {
  		## no need to lookup


  		logp.position <- function(x) {

  		x[,1:2] <- transf(x[,1:2], inv = TRUE)
		## Calculate elevations - note we must expand locations per

		## segment to locations per observation

		elev <- elevation(x[is,1],x[is,2],sun)
		ok <- elev >= ekstrom[1] & elev <= ekstrom[2]
		light.sigma[!ok] <- ekstrom[3]
		## Expected (log scale) light level given tag calibration and the
		## attenuation offset
		att <- calibration(elev)+x[is,3]

		## Compare observed (log) light level to the attenuated expected
		## (log) light level to determine the contribution of each
		## observation to log posterior.
		logp <- dnorm(light, att, light.sigma, log=TRUE)

		## Sum over segments + prior
		sapply(split(logp,segments), sum) + dnorm(x[,3], 0, k.sigma, log = TRUE)
	  }


  	}
  	if (ekstrom.limit == "light") {
  		ekstrom[1:2] <- range(calibration(ekstrom[1:2]))


  		  logp.position <- function(x) {
  		  	x[,1:2] <- transf(x[,1:2], inv = TRUE)
			## Calculate elevations - note we must expand locations per

			## segment to locations per observation

			elev <- elevation(x[is,1],x[is,2],sun)

			## Expected (log scale) light level given tag calibration and the
			## attenuation offset
			att <- calibration(elev)+x[is,3]
			ok <- att >= ekstrom[1] & att <= ekstrom[2]
			light.sigma[!ok] <- ekstrom[3]
			## Compare observed (log) light level to the attenuated expected
			## (log) light level to determine the contribution of each
			## observation to log posterior.
			logp <- dnorm(light, att, light.sigma, log = TRUE)

			## Sum over segments + prior
			sapply(split(logp,segments),sum) + dnorm(x[,3], 0, k.sigma, log=TRUE) } }




  } else {
	  logp.position <- function(x) {

		x[,1:2] <- transf(x[,1:2], inv = TRUE)
	    ## Calculate elevations - note we must expand locations per
	    ## segment to locations per observation
	    elev <- elevation(x[is,1],x[is,2],sun)

	    ## Expected (log scale) light level given tag calibration and the
	    ## attenuation offset
	    att <- calibration(elev)+x[is,3]

	    ## Compare observed (log) light level to the attenuated expected
	    ## (log) light level to determine the contribution of each
	    ## observation to log posterior.
	    logp <- dnorm(light,att,light.sigma,log=TRUE)

	    ## Sum over segments + prior
	    sapply(split(logp,segments),sum)+dnorm(x[,3],0,k.sigma,log=TRUE)
	  }
}
  ##
  ## Behavioural contribution to the log posterior
  ##   - if we specify "distance", then the time difference is assumed 1
  dt <- NULL
  ## Times between observations
  if (behav == "distance")
    dt <- rep(1, nlevels(segments))

  if (behav == "speed") {
    dt <- diff(tapply(as.numeric(day),segments,mean))/dt.scale
    if (any(dt <= 0)) stop("There is a problem with segments, check ", which(dt <= 0))
  }

  if (is.null(dt)) stop("No behavioural specification, need behav as \"distance\" or \"speed\"")

  if(behav.dist=="gamma") {
    alpha <- behav.mean^2/behav.sd^2
    beta <- behav.mean/behav.sd^2
    logp.behavioural <- function(k,x1,z,x2) {
      ## Average speed from x1 to z to x2
      spd <- pmax.int(dist(x1,z)+dist(z,x2), 1e-06)/dt[k]
      dgamma(spd,alpha,beta,log=TRUE)
    }
  }
  else
    {
    log.sigma <- sqrt(log(1+behav.mean^2/behav.sd^2))
    log.mu <- log(behav.mean)-log.sigma^2/2
    logp.behavioural <- function(k,x1,z,x2) {
      ## Average speed from x1 to z to x2
      spd <- (dist(x1,z)+dist(z,x2))/dt[k]
      dnorm(spd,log.mu,log.sigma,log=TRUE)
    }
  }

  ##
  ## Locations to be held fixed
  ##
  n <- nlevels(segments)
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
       ## Light levels prediction
       light.predict=light.predict,
       ## Positional contribution to the log posterior
       logp.position=logp.position,
       ## Behavioural contribution to the log posterior
       logp.behavioural=logp.behavioural,
       ## Locations to be held fixed
       fixed.x=fixed.x,
       proposal.x = proposal.x, proposal.z = proposal.z,
       tm = tapply(day, segments, mean)
       ## Suggested starting points
       #start.x=start.x,
       #start.z=start.z)
       )
}

