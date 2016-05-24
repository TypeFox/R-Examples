brackets <-
function(x1, y1, x2, y2, h=NULL, ticks=0.5, curvature=0.5, type=1, col=1, lwd=1, lty=1, xpd=FALSE){
# Eingabetest
if(!is.numeric(curvature)) stop('curvature must be numeric')
if(!is.numeric(type))      stop('type must be a integer, 1 to 5')
if(!is.logical(xpd))       stop('xpd must be TRUE or FALSE')
if(curvature<0) curvature<- 0
if(curvature>1) curvature<- 1

if(length(ticks)==1) if(is.na(ticks)) ticks<- NULL
if(!is.numeric(ticks) & !is.null(ticks))   stop('ticks must be numeric or NULL')
if(length(ticks)>1){
if(any(duplicated(abs(ticks)))) stop('duplicated ticks')
}
################################################################################
# Dicke berechnen
xm <- mean(c(x1, x2))
ym <- mean(c(y1, y2))
coords <-par('usr')
xr <- abs(diff(coords[1:2]))
yr <- abs(diff(coords[3:4]))
if(is.null(h))  h  <- sqrt(xr*yr)/20
dinps  <- par('pin')
xtox   <- dinps[2]/dinps[1]
rato   <- xr/yr
mysig  <- sign(h)

vx = (x2-x1)/rato/xtox
vy = (y2-y1)*rato*xtox
len = sqrt(vx^2 + vy^2)
ux = vy/len
uy = vx/len
x3 = xm + h * ux
y3 = ym - h * uy
ret <- a_get_relative(c(xm, x3), c(ym, y3))
h   <- mysig*dist(cbind(ret[,1], ret[,2]))

################################################################################
rvalues <- a_get_relative(c(x1, x2), c(y1, y2))
x1 <- rvalues[1,1]
x2 <- rvalues[2,1]
y1 <- rvalues[1,2]
y2 <- rvalues[2,2]

brackets<- a_cb_brackets(phi=curvature, ticks=ticks, type=type)

x <- brackets[1,]*dist(cbind(c(x1, x2), c(y1, y2)))
y <- brackets[2,]*h
rout <- a_rotate(cbind(x,y), atan2(y2-y1, x2-x1))

xout <- rout[,1]+x1
yout <- rout[,2]+y1

native <- a_get_native(xout, yout)
oldxpd <- par('xpd')
par(xpd=xpd)
lines(native[,1], native[,2], type='l', col=col, lwd=lwd, lty=lty, xpd=xpd)
par(xpd=oldxpd)
}
