#
# Ordinary Kriging
#

# Entry points for compiled objects
matrixmult <- function(x, y) {
  n <- length(x[,1])
  p <- length(x[1,])
  m <- length(y[1,])

  return(matrix(.C("mmult", PACKAGE="kriging",
      as.double(c(x)),
      as.double(c(y)),
      as.integer(n),
      as.integer(p),
      as.integer(m),
      z = as.double(rep(0, n*m)))$z, n, m))
}

onedim <- function(x, n) {
  .C("onedimdist", PACKAGE="kriging",
      as.double(x), 
      as.integer(n),
      v = as.double(rep(0, n*n)))$v
}

twodim <- function(x, y, n) {
  .C("twodimdist", PACKAGE="kriging", 
      as.double(x), 
      as.double(y), 
      as.integer(n),
      d = as.double(rep(0, n*n)))$d
}

krig.fit <- function(D, nugget, range, sill, model, n) {
  if(model=="spherical") modelID <- 0
  if(model=="exponential") modelID <- 1
  if(model=="gaussian") modelID <- 2
  return(.C("krigfit", PACKAGE="kriging",
    as.double(D),
    as.double(nugget),
    as.double(range),
    as.double(sill),
    as.integer(modelID),
    as.integer(n),
    a = as.double(rep(1, (n+1)*(n+1))))$a)
}

krig.grid <- function(blx, bly, trx, try, pixels) {
  pixel <- max((trx-blx), (try-bly)) / (pixels-1)
  xpixels <- ceiling((trx-blx)/pixel)
  ypixels <- ceiling((try-bly)/pixel)

  G.grid <- .C("kriggrid", PACKAGE="kriging",
      as.double(blx),
      as.double(bly),
      as.double(pixel),
      as.integer(xpixels),
      as.integer(ypixels),
      gx = as.double(rep(0, xpixels * ypixels)),
      gy = as.double(rep(0, xpixels * ypixels)))

  return(data.frame(x = G.grid$gx, y = G.grid$gy, pixel=pixel))
}

krig.polygons <- function(X, Y, polygons) {
  n <- length(X)
  nlist <- length(polygons)
  npoly <- rep(0, nlist+1)

  polygonsx <- polygonsy <- {}
  for(i in 1:nlist) {
    polygonsx <- c(polygonsx, polygons[[i]][,1])
    polygonsy <- c(polygonsy, polygons[[i]][,2])
    npoly[i+1] <- length(polygonsx)
  }

  G <- .C("krigpolygons", 
    as.integer(n),
    as.double(X),
    as.double(Y),
    as.integer(nlist),
    as.integer(npoly),
    as.double(polygonsx),
    as.double(polygonsy),
    Gn = as.integer(0),
    Gx = as.double(rep(0, n)),
    Gy = as.double(rep(0, n)))

  return(data.frame(x = G$Gx[1:G$Gn], y = G$Gy[1:G$Gn]))
}

krig.pred <- function(x, y, response, Gx, Gy, invA, nugget, range, sill, model, n) {
  if(model=="spherical") modelID <- 0
  if(model=="exponential") modelID <- 1
  if(model=="gaussian") modelID <- 2
  return(.C("krigpred", PACKAGE="kriging",
    as.double(x),
    as.double(y),
    as.double(response),
    as.double(Gx),
    as.double(Gy),
    as.double(invA),
    as.double(nugget),
    as.double(range),
    as.double(sill),
    as.integer(modelID),
    as.integer(n),
    as.integer(length(Gx)),
    as.integer(n+1),
    as.integer(1),
    G.pred = as.double(rep(0, length(Gx))))$G.pred)
}

imageloop <- function(x, y, z, a, b) {
  n <- length(z)
  na <- length(a)
  nb <- length(b)

  loop <- .C("krigimage",
    as.integer(x),
    as.integer(y),
    as.double(z),
    as.integer(n),
    as.integer(a),
    as.integer(b),
    as.integer(na),
    as.integer(nb),
    i = as.double(rep(0, na*nb)),
    j = as.integer(rep(0, na*nb)))

  loop$i[loop$j==0] <- NA
  return(matrix(loop$i, na, nb, byrow=T))
}
