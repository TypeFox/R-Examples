"mkElevationSeg" <-
function(segments,day) {
  ## day - times as POSIXct
  ## segments - indexs separate segments of light curve 
  
  ## Extract components of time (GMT)
  tm <- as.POSIXlt(day,tz="GMT")           
  hh <- tm$hour
  mm <- tm$min
  ss <- tm$sec


  ## Time as Julian day
  jday <- julday(tm)+(hh+(mm+ss/60)/60)/24
  
  ## Time as Julian century
  t <- (jday-2451545)/36525

  ## Geometric mean anomaly for the sun (degrees)
  M <- 357.52911+t*(35999.05029-0.0001537*t)

  ## Equation of centre for the sun (degrees)
  eqcent <- sin(pi/180*M)*(1.914602-t*(0.004817+0.000014*t))+
    sin(pi/180*2*M)*(0.019993-0.000101*t)+
      sin(pi/180*3*M)*0.000289
    
  ## The geometric mean sun longitude (degrees)
  L0 <- 280.46646+t*(36000.76983+0.0003032*t)
  ## Limit to [0,360)
  L0 <- L0%%360 

  ## The true longitude of the sun (degrees)
  lambda0 <- L0 + eqcent
  
  ## The apparent longitude of the sun (degrees)
  omega <- 125.04-1934.136*t
  lambda <- lambda0-0.00569-0.00478*sin(pi/180*omega)
  

  ## The mean obliquity of the ecliptic (degrees)
  seconds <- 21.448-t*(46.815+t*(0.00059-t*(0.001813)))
  obliq0 <- 23+(26+(seconds/60))/60 

  ## The corrected obliquity of the ecliptic (degrees)
  omega <- 125.04-1934.136*t
  obliq <- obliq0 + 0.00256*cos(pi/180*omega)

  ## The eccentricity of earth's orbit
  e <- 0.016708634-t*(0.000042037+0.0000001267*t)
  
  ## The equation of time (minutes of time)
  y <- tan(pi/180*obliq/2)^2
  eqtime <- 180/pi*4*(y*sin(pi/180*2*L0) -
                      2*e*sin(pi/180*M) +
                      4*e*y*sin(pi/180*M)*cos(pi/180*2*L0) -
                      0.5*y^2*sin(pi/180*4*L0) -
                      1.25*e^2*sin(pi/180*2*M))

  ## The sun's declination (radians)
  solarDec <- asin(sin(pi/180*obliq)*sin(pi/180*lambda))
  sinSolarDec <- sin(solarDec)
  cosSolarDec <- cos(solarDec)

  ## ALT
  ## sinSolarDec <- sin(pi/180*obliq)*sin(pi/180*lambda)
  ## cosSolarDec <- cos(asin(sinSolarDec))

  ## Solar time unadjusted for longitude (degrees)
  solarTime <- (hh*60+mm+ss/60+eqtime)/4

  ## Split by segment
  solarTime <- split(solarTime,segments)
  sinSolarDec <- split(sinSolarDec,segments)
  cosSolarDec <- split(cosSolarDec,segments)
  
  function(segment,lon, lat) {
  
    ## Suns hour angle (degrees)
 ## change subtraction to addition to work with -180<->180 convention MDS2Jul03
   
    hourAngle <- solarTime[[segment]]+lon-180
    
    ## Cosine of sun's zenith 
    cosZenith <- sin(pi/180*lat)*sinSolarDec[[segment]]+ 
      cos(pi/180*lat)*cosSolarDec[[segment]]*cos(pi/180*hourAngle)
    
    ## Limit to [-1,1] 
    cosZenith[cosZenith > 1] <- 1
    cosZenith[cosZenith < -1] <- -1
    
    ## Ignore refraction correction
    90-180/pi*acos(cosZenith) 
  }
}



"mkNLPosterior" <-
function(segments,day,light,calib) {

  segments <- unclass(factor(segments))
  
  ## Construct elevation function
  elevation <- mkElevationSeg(segments,day)

  ## Split light levels by segments
  light <- split(light,segments)
  
  ## We return a function that computes the negative log posterior
  function(seg,ps,prv,nxt) {
    ## seg - the segment of the light curve
    ## ps - params for this segment
    ## prv,nxt - params for previous and next segments

    ## Decompose params into (lat, lon) and offset
    #lat <- ps[1]
    #lon <- ps[2]
    lon <- ps[1]
    lat <- ps[2]
    k <- ps[3]

##-  Problem sometimes with elevation - "trying to subscript too many elements"
##- (might be the split)
    ## Compute elevations for this (lat, lon) and segment
    elev <- elevation(seg, lon, lat)

    ## The expected light levels
    lgt <- calib(elev)
    
    ## The log likelihood.
    sigma <- 7
    shape <- (70/30)^2
    rate <- (70/30^2)*10
    eps <- 1.0E-6
    loglik <- sum(dnorm(light[[seg]], k+lgt, 7, log=TRUE),
                  ## Must not allow zero distances for gamma pdf.
                  dgamma(max(eps, old.dist.gc(prv[1:2],ps[1:2])), shape, rate,log=TRUE),
                  dgamma(max(eps, old.dist.gc(ps[1:2],nxt[1:2])), shape, rate,log=TRUE))

    ## Return negative log posterior
    -(log(k.prior(seg, ps))+loglik)
  }
}



"old.dist.gc" <-
function (x1, x2 = NULL) 
{
  if (is.null(x2)) {
    if (!is.matrix(x1) && !is.data.frame(x1)) stop("No points to calculate")
    pt1 <- x1[-1,]
    pt2 <- x1[-nrow(x1),]
    x1 <- pt1
    x2 <- pt2
  }
  pt1 <-  matrix(as.vector(as.matrix(x1)), ncol = 2)
  pt2 <-  matrix(as.vector(as.matrix(x2)), ncol = 2)
    a <- cos(pi/180 * pt1[, 2]) * cos(pi/180 * pt2[, 
        2]) * cos(pi/180 * (pt2[, 1] - pt1[, 1])) + sin(pi/180 * 
        pt1[, 2]) * sin(pi/180 * pt2[, 2])

        6378.137 * acos(pmin(a,1))
}

"old.find.init" <-
function (mask, nseg, nlpost, pars = c("Lon", "Lat", "k")) 
{
    inits <- matrix(0, nseg, length(pars), dimnames = list(NULL, 
        pars))
    xx <- mask$x
    yy <- mask$y
    if (length(xx) != dim(mask$z)[1]) {
        xx <- xx[-1] - diff(xx[1:2])/2
        yy <- yy[-1] - diff(xx[1:2])/2
    }
    xy <- as.matrix(expand.grid(x = xx, y = yy))
   pts <- xy[as.vector(mask$z), ]
    print(paste("k of", nseg))
    for (k in 1:nrow(inits)) {
               
        p <- c(pts[1, 1:2], 0)
        nlp.mode <- nlpost(k, p, p, p)
        p.mode <- p
        for (i in 2:nrow(pts)) {
            p <- c(pts[i, ], 0)
            nlp <- nlpost(k, p, p, p)
            if (nlp < nlp.mode) {
                nlp.mode <- nlp
                p.mode <- p
            }
        }
        cat(k, "\n")
        inits[k, ] <- p.mode
    }
    inits
}

"old.metropolis" <-
function(nlpost,lookup,p0,cov0,start,end, iter=1000,step=100) {
# nlpost<-nlogpost;p0<-ps;cov0<-cov.ps;mat<-summ;iter=20;step=10
  ## Initialize chain - internally we always work with transpose p.
  p <- t(p0)
  
  ## m parameters per twilight and n twilights
  m <- nrow(p)
  n <- ncol(p)

  ## Precalculate the Cholesky decomposition the covariance of the
  ## proposal distribution for each block.
  U <- cov0
  for(j in 1:n) U[,,j] <- chol(U[,,j],pivot = TRUE)## I think you need pivot for 1.7.0 

  ## Record chain as a matrix
  chain.p <- matrix(0,iter,m*n)
  colnames(chain.p) <- paste(colnames(p0),rep(1:n,each=m),sep="")
  
  ## Posterior for each block at initial locations
  l <- rep(0.0, n)
  l[1] <- nlpost(1, p[,1], start, p[,2])
  for(j in 2:(n-1)) {
    l[j] <- nlpost(j,p[,j],p[,j-1],p[,j+1])
  }
  l[n] <- nlpost(n,p[,n],p[,n-1],end)
  
  for(i in 1:iter) {
    for(k in 1:step) {
      
      ## First location compares to start
      pj.k <- p[,1] + rnorm(m) %*% U[,,1]

      ## CHECK this use of test.mask


      ##if(test.mask(masks,1,pj.k[1],pj.k[2])) {
      if (lookup(pj.k[1:2], 1)) {
        lj.k <- nlpost(1,pj.k, start, p[,2])
        ## Test candidate against l for this segment
        if(l[1]-lj.k > log(runif(1))) {
          ## accept candidate
          p[,1] <- pj.k
          l[1] <- lj.k
      

        }
      }
      
      ## Middle locations
      for(j in 2:(n-1))  {
        pj.k <- p[,j] + rnorm(m) %*% U[,,j]
        ##if(test.mask(masks,j,pj.k[1],pj.k[2])) {
        if (lookup(pj.k[1:2], j)) {
          lj.k <- nlpost(j, pj.k, p[,j-1], p[,j+1])
          if(l[j]-lj.k > log(runif(1))) {
            ## accept candidate
            p[,j] <- pj.k
            l[j] <- lj.k

          }
        }
      } 
      
      ## Last location compares to end
     #???pj.k <- p[,n] + rnorm(n,n,m) %*% U[,,n]
	 pj.k <- p[,n] + rnorm(m) %*% U[,,n]
      ##if(test.mask(masks,n,pj.k[1],pj.k[2])) {
      if (lookup(pj.k[2:1], n)) {
        lj.k <- nlpost(n,pj.k,p[,n-1],end)
        if(l[n]-lj.k > log(runif(1))) {
          ## accept candidate
          p[,n] <- pj.k
          l[n] <- lj.k

        }
      }
    
    }
    chain.p[i,] <- p
  
    print(i)
  }

  p <- t(p)
 # names(p) <- names(p0)
  list(p=chain.p,last=p)
 
}

"old.mkLookup" <-
function (x, binArray = TRUE) 
{
    if (!is.list(x)) 
        stop("x must be an image list")
    csize <- c(diff(x$x[1:2]), diff(x$y[1:2]))
    dimXY <- dim(x$z)
    function(xy, segment = 1) {
        if (missing(segment) & binArray) 
            stop("binary arrays require segment")
        if (!is.numeric(segment)) 
            stop("segment must be numeric")
        if (is.vector(xy)) 
            xy <- matrix(xy, ncol = 2, byrow = TRUE)
        if (!is.matrix(xy)) 
            stop("xy must be a matrix or vector")
        xs <- xy[, 1]
        ys <- xy[, 2]
        i <- round((xs - x$x[1] + csize[1]/2)/csize[1])
        j <- round((ys - x$y[1] + csize[2]/2)/csize[2])
        f <- vector(mode(x$z), length(xs))
        k <- (i > 0 & j > 0 & i <= dimXY[1] & j <= dimXY[2])
        if (any(k)) {
            if (binArray) {
                f[k] <- bits(x$z[i[k], j[k], (segment%/%31) + 
                  1], (segment - 1)%%31)
                f == 1
            }
            else {
		f[k] <- x$z[cbind(i[k], j[k])]
		f
		}
        }
        else FALSE
    }
}

k.prior <- 
function (seg, ps) 
dnorm(ps[3], 0, 10)
