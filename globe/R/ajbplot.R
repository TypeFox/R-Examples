
.placedata <- list(
  nedlands=list(lon=115.82,lat=-31.99),
  curtin=list(lon=115.9240, lat=-32.0010),
  perth=list(lon=115.86,lat=-31.95),
  northpole=list(lon=0,lat=90),
  southpole=list(lon=0,lat=-90),
  casey=list(lon=110.524, lat=-66.282),
  mawson=list(lon=62.8667, lat=-67.6000),
  madrid=list(lon=4+25/60, lat=40+17/60),
  aarhus=list(lon=10.21, lat=56.16),
  aalborg=list(lon=9.92, lat=57.05),
  newyorkcity=
  list(lon=-(74+0/60+23/3600),lat=40+42/60+51/3600),
  titanic=list(lon=-(50+14/60), lat=41+46/60),
  pyongyang=list(lon=125.7381, lat=39.0194),
  everest=list(lon=86.9253, lat=27.9881),
  kilimanjaro=list(lon=37.3533, lat=3.0758)
)

place <- function(placename) {
   stopifnot(is.character(placename))
   pn <- tolower(placename)
   if(pn %in% names(.placedata))
     return(.placedata[[pn]])
   stop(paste("Unrecognised placename", sQuote(placename),
              "-- available places are",
              paste(sQuote(names(.placedata)), collapse=", ")),
        call.=FALSE)
}


flatearth <- function(projection=c("atlas", "cylindrical"),
                      gdata, runlen, asp=NULL, ..., do.plot=TRUE){
  if(missing(projection)) 
    projection <- "atlas"
  if(missing(gdata)) {
    gdata <- get("earth")$coords
  }
  if(missing(runlen)) {
    runlen <- get("earth")$runlen
  }
  lonlat <- ensurelonlat(gdata)
  x <- lon <- lonlat$lon
  y <- lat <- lonlat$lat
  xlim <- c(-180,180)
  xat <- seq(-120,120,by=60)
  xlabs <- paste(abs(xat), ifelse(xat < 0, "W", ifelse(xat == 0, "", "E")))
  ylim <- c(-90,90)
  yat <- seq(-80,80,by=20)
  ylabs <- paste(abs(yat), ifelse(yat < 0, "S", ifelse(yat == 0, "", "N")))
  switch(projection,
         cylindrical={
           y <- sin(y * pi/180)
           yat <- sin(yat * pi/180)
           ylim <- c(-1,1)
         },
         atlas={ }
         )
  if(do.plot) {
    if(missing(asp)) {
      plot(x, y, type = "n", axes=FALSE, xlim=xlim,ylim=ylim,
           xlab = "Longitude", ylab = "Latitude")
      axis(1, at=xat,labels=xlabs)
      axis(2, at=yat,labels=ylabs, las=2)
      box()
    } else {
      asp = asp * diff(xlim)/diff(ylim)
      plot(x, y, type = "n", axes=FALSE, xlim=xlim,ylim=ylim,
           xlab = "", ylab = "",
           asp=asp)
      axis(1, at=xat,labels=xlabs, pos=ylim[1])
      axis(2, at=yat,labels=ylabs, las=2, pos=xlim[1])
      title(ylab="Latitude")
      text(mean(xat), ylim[1] - diff(ylim)/5, "Longitude", pos=1)
      lines(xlim[c(1,2,2,1,1)],ylim[c(1,1,2,2,1)])
    }
  }
  
  # remove zeroes from run length vector
  runlen <- runlen[runlen!=0]
  # do not draw line between points numbered breaks[i] and breaks[i]+1
  breaks <- cumsum(runlen)
  # which ones to draw
  s <- seq(x)[-breaks]
  # go
  if(do.plot)
    segments(x[s],y[s],x[s+1],y[s+1], ...)
  result <- cbind(x[s], y[s], x[s+1], y[s+1])
  z <- rep(FALSE, length(x))
  z[breaks] <- TRUE
  attr(result, "piece") <- as.integer(factor(cumsum(z)[-breaks]))
  return(invisible(result))
}

flatpoints <- function(loc, projection=c("atlas", "cylindrical"), ...,
                       do.plot=TRUE) {
  if(missing(projection)) projection <- "atlas"
  lonlat <- ensurelonlat(loc)
  x <- lon <- lonlat$lon
  y <- lat <- lonlat$lat
  switch(projection,
         cylindrical={
           y <- sin(y * pi/180)
         },
         atlas={ }
         )
  if(do.plot)
    points(x, y, ...)
  return(invisible(list(x=x, y=y)))
}

cross <- function(a,b) {
  # 3D vector cross product
  c(a[2] * b[3] - a[3] * b[2],
    - a[1] * b[3] + a[3] * b[1],
    a[1] * b[2] - a[2] * b[1])
}

dot <- function(a, b) {
  # inner product
  sum(a * b)
}

orthogproj <- function(eye=place("nedlands"), top=place("northpole"), loc) {
  #orthogonal projection of each row of 'loc' as seen from 'eye'
  # compute orthogonal triple
  eye <- ensure3d(eye)
  top <- ensure3d(top)
  if(is.matrix(loc)) {
    if(ncol(loc) != 3)
      stop("matrix \"loc\" should have 3 columns")
  } else if(is.vector(loc)) {
    if(length(loc) != 3)
      stop("vector \"loc\" should have length 3, or be an n x 3 matrix")
  }
  a <- - eye/sqrt(sum(eye^2))
  b <- top - a * sum(a * top)
  b <- b/sqrt(sum(b^2))
  c <- cross(a, b)
  # rotation matrix
  rot <- cbind(c, b, a)
  loc %*% rot
}

spatialpos <- function(lon,lat) {
  # spatial position (x,y,z) 
  # (radius of earth = 1)
  # (origin = centre of earth)
  if(missing(lat)) {
    x <- lon
    if(is.list(x)) {
      if(!is.null(x$lon) && !is.null(x$lat)) {
        lon <- x$lon
        lat <- x$lat
      } else {
        lon <- x[[1]]
        lat <- x[[2]]
      }
    } else if(is.matrix(x) && ncol(x) == 2) {
      lon <- x[,1]
      lat <- x[,2]
    } else if(length(x) == 2) {
      lon <- x[1]
      lat <- x[2]
    } else
    stop("Don't know how to extract longitude and latitude from these data")
  }
  theta <- lon * pi/180
  phi  <- lat * pi/180
  cbind(
    cos(phi) * cos(theta),
    cos(phi) * sin(theta),
    sin(phi)
  )
}


ensure3d <- function(x) {
  # Ensure 'x' is a single 3D vector
  if(!is.vector(x) || length(x) != 3)
    x <- as.vector(spatialpos(x))
  if(!is.vector(x) || length(x) != 3)
    stop("Can't convert to 3D point")
  return(x)
}

ensurelonlat <- function(x) {
  if(is.list(x)) {
    if(!is.null(x$lon) && !is.null(x$lat)) {
      lon <- x$lon
      lat <- x$lat
    } else {
      lon <- x[[1]]
      lat <- x[[2]]
    }
  } else if(is.matrix(x) && ncol(x) == 2) {
    lon <- x[,1]
    lat <- x[,2]
  } else if(length(x) == 2) {
    lon <- x[1]
    lat <- x[2]
  } else
  stop("Don't know how to extract longitude and latitude from these data")
  return(list(lon=lon, lat=lat))
}

globeearth <- function(gdata, runlen,
                       eye=place("nedlands"), top=place("northpole"),
                       ..., do.plot=TRUE) {

  if(missing(gdata)) {
    gdata <- get("earth")$coords
  }
  if(missing(runlen)) {
    runlen <- get("earth")$runlen
  }
  
  eye <- ensure3d(eye)
  top <- ensure3d(top)
  
  spos <- spatialpos(gdata[,1],gdata[,2])
  mpos <- orthogproj(eye, top, spos)

  if(do.plot)
    plot(c(-1,1), c(-1,1), type = "n", asp=1, axes=FALSE, xlab="", ylab="")

  x <- mpos[,1]
  y <- mpos[,2]
  ok <- (mpos[,3] < 0)

  # remove initial 0
  runlen <- runlen[runlen!=0]

  breaks <- cumsum(runlen)
  ok[breaks] <- FALSE

  s <- seq(x)[ok]

  if(do.plot) {
    segments(x[s],y[s],x[s+1],y[s+1], ...)
    ## draw globe
    a <- seq(0,2 * pi, length=360)
    lines(cos(a),sin(a),lty=2)
  } 
  result <- cbind(x[s], y[s], x[s+1], y[s+1])
  attr(result, "piece") <- as.integer(factor(cumsum(!ok)[ok]))
  return(invisible(result))
}

globepoints <- function(loc, eye=place("nedlands"), top=place("northpole"), ...,
                        do.plot=TRUE) {
  
  eye <- ensure3d(eye)
  top <- ensure3d(top)

  lonlat <- ensurelonlat(loc)
  spos <- spatialpos(lonlat$lon,lonlat$lat)
  
  mpos <- orthogproj(eye, top, spos)
  x <- mpos[,1]
  y <- mpos[,2]
  ok <- (mpos[,3] < 0)

  result <- list(x=x[ok], y=y[ok])
  if(do.plot)
    points(result, ...)
  return(invisible(result))
}

globelines <- function(loc, eye=place("nedlands"), top=place("northpole"), ...,
                       do.plot=TRUE) {
  
  eye <- ensure3d(eye)
  top <- ensure3d(top)

  lonlat <- ensurelonlat(loc)
  spos <- spatialpos(lonlat$lon,lonlat$lat)
  
  mpos <- orthogproj(eye, top, spos)
  x <- mpos[,1]
  y <- mpos[,2]
  ok <- complete.cases(mpos) & (mpos[,3] < 0)
  n <- nrow(mpos)
  x0 <- x[-n]
  x1 <- x[-1]
  y0 <- y[-n]
  y1 <- y[-1]
  ok <- ok[-n] & ok[-1]
  if(do.plot)
    segments(x0[ok], y0[ok], x1[ok], y1[ok], ...)
  result <- cbind(x0[ok], y0[ok], x1[ok], y1[ok])
  attr(result, "piece") <- as.integer(factor(cumsum(!ok)[ok]))
  return(invisible(result))
}

globedrawlong <- function(lon, eye=place("nedlands"),
                          top=place("northpole"), ...) {
  globelines(expand.grid(lat=c(seq(-90,90,by=1), NA), lon=lon)[,2:1],
             eye, top, ...)
}

globedrawlat <- function(lat, eye=place("nedlands"),
                         top=place("northpole"), ...) {
  globelines(expand.grid(lon=c(seq(-180,180,by=1), NA), lat=lat), eye, top, ...)
}

runifsphere <- function(n) {
  lon <- runif(n,min=-180,max=180)
  lat <- (180/pi) * acos(runif(n, min=-1, max=1))
  lat <- (lat %% 180) - 90
  return(data.frame(lon=lon,lat=lat))
}

runifsphere.wrong <- function(n) {
  lon <- runif(n,min=-180,max=180)
  lat <- runif(n, min=-90,max=90)
  return(data.frame(lon=lon,lat=lat))
}

globearrows <- function(loc, eye=place("nedlands"), top=place("northpole"), len=0.3, ..., do.plot=TRUE) {
  eye <- ensure3d(eye)
  top <- ensure3d(top)
  spos <- spatialpos(loc[,1], loc[,2])
  tailpos <- orthogproj(eye, top, spos)
  headpos <- orthogproj(eye, top, (1 + len) * spos)
  ok <- (tailpos[,3] < 0)
  if(do.plot)
    segments(tailpos[ok,1],tailpos[ok,2],headpos[ok,1],headpos[ok,2], ...)
  result <- cbind(tailpos, headpos)[ok, , drop=FALSE]
  attr(result, "piece") <- as.integer(factor(cumsum(ok)[ok]))
  return(invisible(result))
}

