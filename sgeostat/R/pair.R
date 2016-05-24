assign("pair",
function (point.obj,num.lags=10,type='isotropic',theta=0,dtheta=5,maxdist) {
  if (!inherits(point.obj,"point")) stop('Point.Obj must be of class, "point".\n')

  if (length(type) != 1) stop('Length of "type" must be 1.\n')
  if ((type != "isotropic") & (type != 'anisotropic'))
    stop('Type must be "isotropic" or "anisotropic".\n')

  if (missing(maxdist)) maxdist <- -1
  if (type[1] == 'isotropic') o.pair <- pair.iso(point.obj,num.lags,maxdist)
  else o.pair <- pair.aniso(point.obj,num.lags,theta,dtheta,maxdist)

  cat("\n")  
  class(o.pair) <- "pair"
  return(o.pair)
})


#************ pair.iso *******************
# creates pointers for the n choose 2 possible pairs of points...
#******************************************

assign("pair.iso",
function(point.obj,num.lags,maxdist) {
# we only need some of the pairs...
  max.entered <- TRUE
  if (maxdist < 0) {
    max.entered <- FALSE
    maxdist <- ((max(point.obj$x)-min(point.obj$x))^2 + (max(point.obj$y)-min(point.obj$y))^2)^0.5
  }

  from <- numeric(0)
  to <- numeric(0)
  n <- length(point.obj$x)
  if (n > 1000) cat('Creating from, to vectors')
  for (i in 1:(n-1)) {
    cat ('.')
     
    candidates <- (i+1):n
    candidates <- candidates[(point.obj$x[candidates] > point.obj$x[i]-maxdist) &
                            (point.obj$x[candidates] < point.obj$x[i]+maxdist) &
                            (point.obj$y[candidates] > point.obj$y[i]-maxdist) &
                            (point.obj$y[candidates] < point.obj$y[i]+maxdist)]

    from <- c(from,rep(i,length(candidates)))
    to <- c(to,candidates)
  }
  if (n > 1000) cat('\n')

# calculate the distance between all possible pairs...
  if (n > 1000) cat('Calculating distances...')
  dist <- sqrt( (point.obj$x[from]-point.obj$x[to])^2 + 
                (point.obj$y[from]-point.obj$y[to])^2 )
  if (n > 1000) cat('\n')

# apply the maximum distance cutoff, if specified...
  if (maxdist < 0) maxdist <- max(dist,na.rm=TRUE)
  if (max.entered) {
    if (n > 1000) cat('Applying maxdist criterion...')
    from <- from[dist<=maxdist]
    to   <- to[dist<=maxdist]
    dist <- dist[dist<=maxdist]
    if (n > 1000) cat('\n')
  }
  else maxdist <- max(dist)

# create the vector to use to "cut" the bins...
  bins.cut <- seq(0,maxdist,maxdist/num.lags)

# create the vector of bins center points (for plotting)...
#  bins.cent <- list()
# revision 11/4/99 rsb
  bins.cent <- numeric(length(bins.cut)-1)
  for (i in 1:(length(bins.cut)-1))
    bins.cent[i] <- bins.cut[i]+(bins.cut[i+1]-bins.cut[i])/2

# cut the data into lag bins...
  if (n > 1000) cat('Cutting distances into bins...')
#  lags <- cut(dist,bins.cut)
  lags <- cut(dist,bins.cut,labels=c(1:num.lags))
  if (n > 1000) cat('\n')

# Journel says that you should have at least 30 pairs of points
# give a warning if not
  for(i in 1:num.lags) {
    if (length(lags[lags==1]) < 30) 
      cat(paste('NOTE: Number of pairs in lag ',i,': ',length(lags[lags==i]),'\n',collapse=""))
  }

  pair <- list(from=from,to=to,lags=lags,dist=dist,bins=bins.cent)
  attr(pair,"type") <- 'isotropic'
  attr(pair,"theta") <- NULL
  attr(pair,"dtheta") <- NULL

  return(pair)
})

#
#************ pair.aniso *******************
# Creates pairs that fall within a given direction and angle
#********************************************
assign("pair.aniso",
function(point.obj,num.lags,theta,dtheta,maxdist) {

# we only need some of the pairs...
  from <- NULL
  to <- NULL
  n <- length(point.obj$x)
  if (n > 1000) cat('Creating from, to vectors')
  for (i in 1:(n-1)) {
    if (i%%5==0) cat ('.')
     
    candidates <- (i+1):n
    candidates <- candidates[point.obj$x[candidates] > point.obj$x[i]-maxdist &
                            point.obj$x[candidates] < point.obj$x[i]+maxdist &
                            point.obj$y[candidates] > point.obj$y[i]-maxdist &
                            point.obj$y[candidates] < point.obj$y[i]+maxdist]

    from <- c(from,rep(i,length(candidates)))
    to <- c(to,candidates)
  }
  cat('\n')
  if (n > 1000) cat('\n')
  
# look both ways...
  xx <- from
  from <- c(from,to)
  to   <- c(to,xx)

# Calculate the distance...
  dist <- sqrt( (point.obj$x[from]-point.obj$x[to])^2 + 
                (point.obj$y[from]-point.obj$y[to])^2 )

# Apply the maximum distance criteria, if entered...
  if (maxdist < 0) maxdist <- max(dist,na.rm=TRUE)
  from <- from[dist<=maxdist]
  to   <- to[dist<=maxdist]
  dist <- dist[dist<=maxdist]

# calc the angle between pairs
  angle <- calcangle(point.obj$x[from],point.obj$y[from],
                 point.obj$x[to],point.obj$y[to])

# if two points have the same location, they will be NA's
  to    <- to[!is.na(angle)]
  from  <- from[!is.na(angle)]
  dist  <- dist[!is.na(angle)]
  angle <- angle[!is.na(angle)]

# convert theta and dtheta to radians...
  theta.rad <- 2*pi*theta/360
  dtheta.rad <- 2*pi*dtheta/360

# Get the angle to look for data...
# need to be careful around angle 0...
  startangle <- theta.rad-dtheta.rad
  endangle <- theta.rad+dtheta.rad
  if(startangle<0) startangle <- 2*pi+startangle
  else if (endangle>2*pi) endangle <- endangle-2*pi

# Apply the angle criteria...
  if (startangle>endangle) {
    from <- from[angle>startangle | angle<endangle]
    to <- to[angle>startangle | angle<endangle]
    dist <- dist[angle>startangle | angle<endangle]
    angle <- angle[angle>startangle | angle<endangle]
  }
  else {
    from <- from[angle>startangle & angle<endangle]
    to <- to[angle>startangle & angle<endangle]
    dist <- dist[angle>startangle & angle<endangle]
    angle <- angle[angle>startangle & angle<endangle]
  }

# create the vector to use to "cut" the bins...
  bins.cut <- seq(0,max(dist,na.rm=TRUE),max(dist,na.rm=TRUE)/num.lags)

# create the vector of bin center points (for plotting)...
  bins.cent <- NULL
  for (i in 1:(length(bins.cut)-1))
    bins.cent[i] <- bins.cut[i]+(bins.cut[i+1]-bins.cut[i])/2

# cut the pairs into lags...
  lags <- cut(dist,bins.cut,labels=c(1:num.lags))

# ?? says that there should be at least 30 pairs in each bin.
# Give a warning otherwise...
  for(i in 1:num.lags) {
    if (length(lags[lags==i]) < 30) 
      cat(paste('NOTE: Number of pairs in lag ',i,': ',length(lags[lags==i]),'\n',collapse=""))
  }

  pair <- list(from=from,to=to,lags=lags,dist=dist,bins=bins.cent)
  attr(pair,"type") <- 'anisotropic'
  attr(pair,"theta") <- format(theta)
  attr(pair,"dtheta") <- format(dtheta)

  return(pair)
})


#**************
assign("calcangle",
function(x1,y1,x2,y2) {

  xdoty <- x2-x1
  lenvect <- sqrt( (x2-x1)^2 + (y2-y1)^2 )
  angle <- acos( xdoty / lenvect )

  angle <- ifelse( y2 < y1, (2*pi)-angle, angle)

  return(angle)
})

#
#
#************ pair.newangle *******************
# Takes an anisotopic pair object and makes a
# isotropic version.  Don't have to search the thing.
#********************************************
assign("pair.newangle",
function(point.obj, pair.obj,num.lags=10,theta=0,dtheta=5,maxdist) {

  if (!inherits(point.obj,"point")) stop('Point.obj must be of class, "point".\n')
  if (!inherits(pair.obj,"pair")) stop('Pair.obj must be of class, "pair".\n')
  if (attr(pair.obj,"type") != 'isotropic') stop('Pair.obj must be isotropic.\n')

  from <- pair.obj$from
  to   <- pair.obj$to

# look both ways...
  xx <- from
  from <- c(from,to)
  to   <- c(to,xx)

# Calculate the distance...
  dist <- sqrt( (point.obj$x[from]-point.obj$x[to])^2 + 
                (point.obj$y[from]-point.obj$y[to])^2 )

# Apply the maximum distance criteria, if entered...
  if (maxdist < 0) maxdist <- max(dist,na.rm=TRUE)
  from <- from[dist<=maxdist]
  to   <- to[dist<=maxdist]
  dist <- dist[dist<=maxdist]

# calc the angle between pairs
  angle <- calcangle(point.obj$x[from],point.obj$y[from],
                 point.obj$x[to],point.obj$y[to])

# if two points have the same location, they will be NA's
  to    <- to[!is.na(angle)]
  from  <- from[!is.na(angle)]
  dist  <- dist[!is.na(angle)]
  angle <- angle[!is.na(angle)]

# convert theta and dtheta to radians...
  theta.rad <- 2*pi*theta/360
  dtheta.rad <- 2*pi*dtheta/360

# Get the angle to look for data...
# need to be careful around angle 0...
  startangle <- theta.rad-dtheta.rad
  endangle <- theta.rad+dtheta.rad
  if(startangle<0) startangle <- 2*pi+startangle
  else if (endangle>2*pi) endangle <- endangle-2*pi

# Apply the angle criteria...
  if (startangle>endangle) {
    from <- from[angle>startangle | angle<endangle]
    to <- to[angle>startangle | angle<endangle]
    dist <- dist[angle>startangle | angle<endangle]
    angle <- angle[angle>startangle | angle<endangle]
  }
  else {
    from <- from[angle>startangle & angle<endangle]
    to <- to[angle>startangle & angle<endangle]
    dist <- dist[angle>startangle & angle<endangle]
    angle <- angle[angle>startangle & angle<endangle]
  }

# create the vector to use to "cut" the bins...
  bins.cut <- seq(0,max(dist,na.rm=TRUE),max(dist,na.rm=TRUE)/num.lags)

# create the vector of bin center points (for plotting)...
  bins.cent <- NULL
  for (i in 1:(length(bins.cut)-1))
    bins.cent[i] <- bins.cut[i]+(bins.cut[i+1]-bins.cut[i])/2

# cut the pairs into lags...
  lags <- cut(dist,bins.cut)

# ?? says that there should be at least 30 pairs in each bin.
# Give a warning otherwise...
  for(i in 1:num.lags) {
    if (length(lags[lags==i]) < 30) 
      cat(paste('NOTE: Number of pairs in lag ',i,': ',length(lags[lags==i]),'\n',collapse=""))
  }

  pair <- list(from=from,to=to,lags=lags,dist=dist,bins=bins.cent)
  attr(pair,"type") <- 'anisotropic'
  attr(pair,"theta") <- format(theta)
  attr(pair,"dtheta") <- format(dtheta)

  attr(pair,"class") <- "pair"
  return(pair)

})
