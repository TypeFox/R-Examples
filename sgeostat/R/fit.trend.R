
assign("fit.trend",
function (point.obj,at,np=2,plot.it=TRUE) {
  if (!inherits(point.obj,"point")) stop('Point.Obj must be of class, "point".\n')

  if(missing(at)) stop('Must enter at least one attribute.\n')

  z <- point.obj[[match(at,names(point.obj))]]

  X <- point.obj$x
  Y <- point.obj$y
  P <- ((np + 1) * (np + 2))/2
  mat <- matrix(1, length(X), P)
  dimnames(mat) <- list(NULL, c("const", rep("", P - 1)))
  ip <- 0
  if(np > 0)
    for(i in (0:np)) {
      for(j in (0:(np - i))) {
        ip <- ip + 1
        mat[, ip] <- X^j * Y^i
        dimnames(mat)[[2]][ip] <- paste("x^", j, " y^", 
        i, sep = "")
      }
    }
  if(length(z) != length(X))
    stop("lengths of x and z must match")
  z.qr <- qr(mat)
  R <- qr.R(z.qr)
  beta <- qr.coef(z.qr, z)
  W <- qr.resid(z.qr, z)
  ts <- structure(list(beta = beta, R = R,  np = np,
    x = X, y = Y, z = z, residuals = W), class = "trend.surface")

  if(plot.it) {
    n <- 25 # a density of 25x25
    tr.mat<-trend.matrix(ts,min(X),max(X),min(Y),max(Y),30)
#    points(perspp(X,Y,z,
#      persp(trend.matrix(ts,min(X),max(X),min(Y),max(Y),30))))
    contour(tr.mat$x,tr.mat$y,tr.mat$z)
    points(X,Y)
  }
    
  return(ts)
})

assign("trend.matrix",
function (ts.obj,xl,xu,yl,yu,n) {

# Adapted from Ripley's Trmat function...
  if(!inherits(ts.obj, "trend.surface"))
    stop("object not a fitted trend surface")
  dx <- (xu - xl)/n
  dy <- (yu - yl)/n
  x <- seq(xl, xu, dx)
  y <- seq(yl, yu, dx)
  z <- matrix(nrow = length(x), ncol = length(y))
  for(i in 1:length(y))
    z[, i] <- trend.value(ts.obj, x, rep(y[i], length(x)))
  invisible(list(x = x, y = y, z = z))
})

assign("trend.value",
function (ts.obj,x,y) {

# Adapted from Ripley's Trval function...
  if(length(x) != length(y))
    stop("lengths of x and y must match")
  degree <- ts.obj$np
  P <- ((degree + 1) * (degree + 2))/2
  mat <- matrix(1, length(x), P)
  ip <- 0
  if(degree > 0)
    for(i in (0:degree)) {
      for(j in (0:(degree - i))) {
        ip <- ip + 1
        mat[, ip] <- x^j * y^i
      }
    }
  as.vector(mat %*% ts.obj$beta)

})
# Versuch ...
"prediction.matrix"<-function (point.obj, v1, var.mod.object,
			       xl,xu,yl,yu,n,maxdist=NULL){
  if(!inherits(point.obj, "point"))
    stop("object not a point object")
  if(!inherits(var.mod.object, "variogram.model"))
    stop("object not a variogram.model object")
  dx <- (xu - xl)/n
  dy <- (yu - yl)/n
  x <- seq(xl, xu, dx)
  y <- seq(yl, yu, dx)
  z <- matrix(nrow = length(x), ncol = length(y))
  
  for(i in 1:length(y)){
    for(j in 1:length(x)){
      predpnt<-point(data.frame(x=x[j],y=y[i]))
      predpnt<- krige(predpnt,point.obj,v1,var.mod.object,maxdist)
      z[j, i]<-predpnt$zhat
    }
    invisible(list(x = x, y = y, z = z))
  }
}
