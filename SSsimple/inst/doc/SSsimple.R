### R code from vignette source 'SSsimple.Snw'

###################################################
### code chunk number 1: SSsimple.Snw:122-123
###################################################
options(width=75)


###################################################
### code chunk number 2: a1
###################################################
library(SSsimple)
F <- 1
H <- matrix(1) ### a 1 x 1 matrix, implying to SS.sim() that n=1 and d=1
Q <- 1
R <- 0
tt <- 1000
set.seed(999)
xss <- SS.sim(F=F, Q=Q, H=H, R=R, length.out=tt, beta0=0)


###################################################
### code chunk number 3: SSsimple.Snw:144-145
###################################################
plot( xss$Z, type="l", main="Random Walk, Sim" )


###################################################
### code chunk number 4: a2
###################################################
F <- matrix( c( 0.65, 0.3, 0.3, 0.65 ), 2, 2 ) ;
H <- matrix( c(0.4, 0.3), 1, 2 ) ### a 1 x 2 matrix, so SS.sim() knows n=1 and d=2 ;
Q <- 0.2 ;
R <- 5 ;
tt <- 180 ;
set.seed(999) ;
xss <- SS.sim(F=F, Q=Q, H=H, R=R, length.out=tt, beta0=0) ;


###################################################
### code chunk number 5: SSsimple.Snw:168-169
###################################################
plot( xss$Z, type="l", main="SPM minus FPM, Sim" )


###################################################
### code chunk number 6: a3
###################################################
P0 <- diag(Q, 2) %*% solve(diag(1, 2) - t(F) %*% F)
xslv <- SS.solve(Z=xss$Z, F=F, Q=Q, H=H, R=R, length.out=tt, P0=P0, beta0=0)
Z.hat <- t(H %*% t(xslv$B.apri))
sqrt( mean( (xss$Z - Z.hat)^2 ) )


###################################################
### code chunk number 7: SSsimple.Snw:185-186
###################################################
options(width=55)


###################################################
### code chunk number 8: SSsimple.Snw:190-199
###################################################
par( mfrow=c(1,2) )
plot(xss$Beta[ , 1], type="l", ylim=range(xss$Beta), col="red", 
   ylab="Humor Index (True is Heavy Line)", 
   main="Humor: True State and Posterior Est State", lwd=4)
points(xslv$B.apos[ , 1], type="l", ylim=range(xslv$B.apos), col="red")
plot(xss$Beta[ , 2], type="l", ylim=range(xss$Beta), col="blue", 
   ylab="Happiness Index (True is Heavy Line)", 
   main="Happiness: True State and Posterior Est State", lwd=4)
points(xslv$B.apos[ , 2], type="l", ylim=range(xslv$B.apos), col="blue")


###################################################
### code chunk number 9: a4
###################################################
xid <- SS.ID( xss$Z, d=2 )
xslv <- SS.solve(Z=xss$Z, F=xid$F, Q=xid$Q, H=xid$H, R=xid$R, 
   length.out=tt, P0=P0, beta0=0)
Z.hat.2 <- t(xid$H %*% t(xslv$B.apri))
sqrt( mean( (xss$Z - Z.hat.2)^2 ) )


###################################################
### code chunk number 10: SSsimple.Snw:216-221
###################################################
par( mfrow=c(1,2) )
plot( xss$Z, type="l", main="SPM - FPM, Prior Est Gold, True Hypers", lwd=3 )
points( Z.hat, type="l", lwd=3, col="gold" )
plot( xss$Z, type="l", main="SPM - FPM, Prior Est Gold, IDed Hypers", lwd=3 )
points( Z.hat.2, type="l", lwd=3, col="gold" )


###################################################
### code chunk number 11: a5
###################################################
d <- 4
n <- 1
F <- matrix(0, d, d)
F[ 1:2, 1:2 ] <- c(0.65, 0.3, 0.3, 0.65)
F[ 3, 3 ] <- 1
F[ 4, 4 ] <- 0.99
eigen(F) #### check stability
tt <- 180

H.tv <- list()
for(i in 1:tt) {
   H.tv[[i]] <- matrix( c(0.4, 0.3, 0, 1), n, d )
   if( i >= 60 & i < 90 ) { H.tv[[i]][ , 3] <- 9 }
   if( i >= 120 & i < 150 ) { H.tv[[i]][ , 3] <- -9 }
}
	
Q <- diag( 0.2, d )
Q[ 3, 3 ] <- 0
Q[ 4, 4 ] <- 1/100
R <- 5

beta0 <- c(0, 0, 1, 0)

set.seed(999)
xss <- SS.sim.tv(F=F, Q=Q, H=H.tv, R=R, length.out=tt, beta0=beta0)


###################################################
### code chunk number 12: SSsimple.Snw:261-266
###################################################
plot( xss$Z, type="l", main="SPM minus FPM, with Video, Sim" ) ;
abline(v=60, col="orange", lwd=2) ;
abline(v=91, col="orange", lwd=2) ;
abline(v=120, col="blue", lwd=2) ;
abline(v=151, col="blue", lwd=2)


###################################################
### code chunk number 13: a6
###################################################
x <- I(0:20) / 20
n <- length(x)
H <- H.omega.sincos(x, c(1,2))
F <- 0.999 ; Q <- 0.1
## the following line constructs the matrix of Euclidean distances btwn locations
D <- abs( tcrossprod(x, rep(1, n)) - tcrossprod(rep(1, n), x) )
R <- exp(-3 * D)

set.seed(999)
xss <- SS.sim(F=F, Q=Q, H=H, R=R, length.out=100, beta0=0)

xdom <- I(0:100) / 100
Hdom <- H.omega.sincos(xdom, c(1,2))


###################################################
### code chunk number 14: a7 (eval = FALSE)
###################################################
## for(i in 1:tt) {
##    plot(x, xss$Z[i, ], ylim=range(xss$Z), main=i)
##    points( xdom, Hdom %*% xss$Beta[i,], type="l" )
##    Sys.sleep(0.1)
## }


###################################################
### code chunk number 15: SSsimple.Snw:307-310
###################################################
i <- 100
   plot(x, xss$Z[i, ], ylim=range(xss$Z), main=i)
   points( xdom, Hdom %*% xss$Beta[i,], type="l" )


###################################################
### code chunk number 16: a10
###################################################
x <- rep( 0:10 / 10, 11 )
y <- rep( 0:10 / 10, each=11 )
n <- length(x)

Hx <- H.omega.sincos( x, c(1,2,3)*pi / 2 )
Hy <- H.omega.sincos( y, c(1,2,3)*pi / 2 )

## Construct the tensor bases expansion over Omega
H <- matrix(NA, nrow(Hx), ncol(Hx)*ncol(Hy))
k <- 0
for(i in 1:ncol(Hx)) {
   for(j in 1:ncol(Hy)) {
      k <- k+1
      H[ , k]  <- Hx[ ,i ] * Hy[ , j]
   }
}

Dx <- tcrossprod(x, rep(1, n)) - tcrossprod(rep(1, n), x)
Dy <- tcrossprod(y, rep(1, n)) - tcrossprod(rep(1, n), y)
D <- sqrt( Dx^2 + Dy^2 ) ; ## the Euclidean distance matrix
R <- exp(-3 * D)

xss <- SS.sim( 0.99, H, 1/2, R,   500, rep(0, ncol(H)) )


###################################################
### code chunk number 17: a11 (eval = FALSE)
###################################################
## for(i in 1:nrow(xss$Z)) {
##    plot(x, y, cex=(xss$Z[i ,]-min(xss$Z))/30, main=i)
##    Sys.sleep(0.1)
## }


###################################################
### code chunk number 18: SSsimple.Snw:365-367
###################################################
i <- 100
plot(x, y, cex=(xss$Z[i ,]-min(xss$Z))/30, main=i) ;


###################################################
### code chunk number 19: a12
###################################################
#################### raster map, interpolation to arb locs
x.grid <- 0:100 / 100
y.grid <- 0:100 / 100
xdom <- rep( x.grid, length(y.grid) )
ydom <- rep( y.grid, each=length(x.grid) )

Hx.dom <- H.omega.sincos( xdom, c(1,2,3)*pi / 2 )
Hy.dom <- H.omega.sincos( ydom, c(1,2,3)*pi / 2 )

Hdom <- matrix(NA, nrow(Hx.dom), ncol(Hx.dom)*ncol(Hy.dom))
k <- 0
for(i in 1:ncol(Hx.dom)) {
   for(j in 1:ncol(Hy.dom)) {
      k <- k+1
      Hdom[ , k]  <- Hx.dom[ ,i ] * Hy.dom[ , j]
   }
}

Dx.dom <- tcrossprod(x, rep(1, length(xdom))) - tcrossprod(rep(1, n), xdom)
Dy.dom <- tcrossprod(y, rep(1, length(ydom))) - tcrossprod(rep(1, n), ydom)
D.dom <- sqrt( Dx.dom^2 + Dy.dom^2 )
R0 <- exp(-3 * D.dom)
bb <- solve(R) %*% R0


###################################################
### code chunk number 20: a13 (eval = FALSE)
###################################################
## for(i in 1:nrow(xss$Z)) {
##    z.hat <- H %*% xss$Beta[i, ]
##    z.0 <- Hdom %*% xss$Beta[i, ]
##    z.tilde <- z.0 + t( t(xss$Z[i, ] - z.hat) %*% bb )
##    Z.mx <- matrix( z.tilde, length(y.grid), length(x.grid) )
##    image(x.grid, y.grid, Z.mx, zlim=range(xss$Z), main=i, col=heat.colors(10000))
##    points(x, y, cex=(xss$Z[i ,]-min(xss$Z))/30)
##    Sys.sleep(0.1)
## }


###################################################
### code chunk number 21: myFigX3
###################################################
png("myFigX3.png")
i <- 100
z.hat <- H %*% xss$Beta[i, ]
z.0 <- Hdom %*% xss$Beta[i, ]
z.tilde <- z.0 + t( t(xss$Z[i, ] - z.hat) %*% bb )
Z.mx <- matrix( z.tilde, length(y.grid), length(x.grid) )
image(x.grid, y.grid, Z.mx, zlim=range(xss$Z), main=i, col=heat.colors(10000))
points(x, y, cex=(xss$Z[i ,]-min(xss$Z))/30)
dev.off()


###################################################
### code chunk number 22: SSsimple.Snw:440-441
###################################################
library(maps)


###################################################
### code chunk number 23: a14
###################################################
## library(maps)
data(SS_O3) #### two dataframes, Z and locs
Z <- SS_O3$Z
locs <- SS_O3$locs
xdate <- row.names(Z)
x <- locs[ ,1]
y <- locs[ ,2]
Z <- as.matrix(Z)
tt <- nrow(Z)
n <- ncol(Z)
Dx <- tcrossprod(x, rep(1, n)) - tcrossprod(rep(1, n), x)
Dy <- tcrossprod(y, rep(1, n)) - tcrossprod(rep(1, n), y)
D <- sqrt( Dx^2 + Dy^2 )


###################################################
### code chunk number 24: a15 (eval = FALSE)
###################################################
## for(i in 1:tt) {
## 	plot(x, y, cex=(Z[i ,]-min(Z))/30, main=xdate[i])
## ##	map("state", "california", add=TRUE)
## 	Sys.sleep(1/5)
## }


###################################################
### code chunk number 25: SSsimple.Snw:468-471
###################################################
i <- 600
plot(x, y, cex=(Z[i ,]-min(Z))/30, main=xdate[i])
map("state", "california", add=TRUE)


###################################################
### code chunk number 26: a16
###################################################
Q <- 1
F <- 1
R <- 1
H <- matrix( 1, n, 1 )

xslv <- SS.solve(Z=Z, F=F, Q=Q, H=H, R=R, length.out=tt, P0=10^5, beta0=0)

Z.hat <- t(H %*% t(xslv$B.apri))
sqrt( mean( ( Z - Z.hat )[10:tt, ]^2 ) )


###################################################
### code chunk number 27: a17
###################################################
ux <- I(1:3)*pi / 20
uy <- I(1:3)*pi / 20
Hx <- H.omega.sincos( x, ux )
Hy <- H.omega.sincos( y, uy )
H <- cbind( rep(1,n), Hx, Hy )
R <- exp( -0.11 * D )
Q <- 1
F <- 1

xslv <- SS.solve(Z=Z, F=F, Q=Q, H=H, R=R, length.out=tt, P0=10^5, beta0=0)

Z.hat <- t(H %*% t(xslv$B.apri))
sqrt( mean( ( Z - Z.hat )[10:tt, ]^2 ) )


###################################################
### code chunk number 28: a18
###################################################
x.grid <- seq( min(x)-0.5, max(x)+0.5, length=100 )
y.grid <- seq( min(y)-0.5, max(y)+0.5, length=100 )
xdom <- rep( x.grid, length(y.grid) )
ydom <- rep( y.grid, each=length(x.grid) )
Hx.dom <- H.omega.sincos( xdom, ux )
Hy.dom <- H.omega.sincos( ydom, uy )
Hdom <- cbind( rep(1, length(xdom)), Hx.dom, Hy.dom )


###################################################
### code chunk number 29: a19 (eval = FALSE)
###################################################
## for(i in 1:nrow(Z)) {
##    Z.mx <- matrix( Hdom %*% xslv$B.apri[i, ], length(y.grid), length(x.grid) )
##    image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000))
##    points(x, y, cex=(Z[i ,]-min(Z))/30)
##    map("state", "california", add=TRUE)
##    Sys.sleep(0.1)
## }


###################################################
### code chunk number 30: myFigX4
###################################################
png("myFigX4.png")
i=600
Z.mx <- matrix( Hdom %*% xslv$B.apri[i, ], length(y.grid), length(x.grid) )
image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000))
points(x, y, cex=(Z[i ,]-min(Z))/30)
map("state", "california", add=TRUE)
dev.off()


###################################################
### code chunk number 31: a20
###################################################
####################### tensor bases
H <- matrix(NA, nrow(Hx), ncol(Hx)*ncol(Hy))
k <- 0
for(i in 1:ncol(Hx)) {
   for(j in 1:ncol(Hy)) {
      k <- k+1
      H[ , k]  <- Hx[ ,i ] * Hy[ , j]
   }
}
H <- cbind( rep(1,n), H ) ### add intercept
R <- exp( -0.11 * D )
Q <- 1
F <- 1

xslv <- SS.solve(Z=Z, F=F, Q=Q, H=H, R=R, length.out=tt, P0=10^5, beta0=0)

Z.hat <- t(H %*% t(xslv$B.apri))
sqrt( mean( ( Z - Z.hat )[10:tt, ]^2 ) )

####### new sites
Hdom <- matrix(NA, nrow(Hx.dom), ncol(Hx.dom)*ncol(Hy.dom))
k <- 0
for(i in 1:ncol(Hx.dom)) {
   for(j in 1:ncol(Hy.dom)) {
      k <- k+1
      Hdom[ , k]  <- Hx.dom[ ,i ] * Hy.dom[ , j]
   }
}
Hdom <- cbind( rep(1, length(xdom)), Hdom )


###################################################
### code chunk number 32: a21 (eval = FALSE)
###################################################
## for(i in 1:nrow(Z)) { ;
##    Z.mx <- matrix( Hdom %*% xslv$B.apri[i, ], length(y.grid), length(x.grid) ) ;
##    image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000)) ;
##    points(x, y, cex=(Z[i ,]-min(Z))/30) ;
##    map("state", "california", add=TRUE) ;
##    Sys.sleep(0.1) ;
## }


###################################################
### code chunk number 33: myFigX2
###################################################
png("myFigX2.png")
i=600
Z.mx <- matrix( Hdom %*% xslv$B.apri[i, ], length(y.grid), length(x.grid) )
image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000))
points(x, y, cex=(Z[i ,]-min(Z))/30)
map("state", "california", add=TRUE)
dev.off()


###################################################
### code chunk number 34: a22
###################################################
H <- diag(1, n)
R <- diag(1, n)
d <- ncol(H)
F <- diag(1, d)
Q <- 1

xslv <- SS.solve(Z=Z, F=F, Q=Q, H=H, R=R, length.out=tt, P0=10^5, beta0=0)

Z.hat <- t(H %*% t(xslv$B.apri))
sqrt( mean( ( Z - Z.hat )[10:tt, ]^2 ) )

x.grid <- seq( min(x)-0.5, max(x)+0.5, length=300 )
y.grid <- seq( min(y)-0.5, max(y)+0.5, length=300 )
xdom <- rep( x.grid, length(y.grid) )
ydom <- rep( y.grid, each=length(x.grid) )

Dx.dom <- tcrossprod(x, rep(1, length(xdom))) - tcrossprod(rep(1, n), xdom)
Dy.dom <- tcrossprod(y, rep(1, length(ydom))) - tcrossprod(rep(1, n), ydom)
D.dom <- t(  sqrt( Dx.dom^2 + Dy.dom^2 )  )
xmin <- apply(D.dom, 1, min)
xmin.mx <- matrix( xmin, nrow(D.dom), ncol(D.dom) )
Hdom <- matrix( as.integer( D.dom == xmin.mx ), nrow(D.dom), ncol(D.dom) )
rm( Dx.dom, Dy.dom, D.dom )


###################################################
### code chunk number 35: a23 (eval = FALSE)
###################################################
## for(i in 1:nrow(Z)) {
##    Z.mx <- matrix( Hdom %*% xslv$B.apri[i, ], length(y.grid), length(x.grid) )
##    image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000))
##    points(x, y, cex=(Z[i ,]-min(Z))/30)
##    map("state", "california", add=TRUE)
##    Sys.sleep(0.1)
## }


###################################################
### code chunk number 36: myFigX1
###################################################
png("myFigX1.png")
i <- 600
Z.mx <- matrix( Hdom %*% xslv$B.apri[i, ], length(y.grid), length(x.grid) )
image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000))
points(x, y, cex=(Z[i ,]-min(Z))/30)
map("state", "california", add=TRUE)
dev.off()


###################################################
### code chunk number 37: a24
###################################################
xid <- SS.ID( Z + rnorm(tt*ncol(Z), 0, 0.1) , d=7, rsN <- c(3, 3, 350) ) ;
xslv <- SS.solve(Z=Z, F=xid$F, Q=xid$Q, H=xid$H, R=xid$R, length.out=tt, 
   P0=10^5, beta0=0) ;
Z.hat <- t(xid$H %*% t(xslv$B.apri)) ;
sqrt( mean( ( Z - Z.hat )[10:tt, ]^2 ) )


###################################################
### code chunk number 38: a25
###################################################
data(SS_O3) #### two dataframes, Z and locs
Z <- SS_O3$Z
locs <- SS_O3$locs
xdate <- row.names(Z)
x <- locs[ ,1]
y <- locs[ ,2]
Z <- as.matrix(Z)
tt <- nrow(Z)
n <- ncol(Z)
Dx <- tcrossprod(x, rep(1, n)) - tcrossprod(rep(1, n), x)
Dy <- tcrossprod(y, rep(1, n)) - tcrossprod(rep(1, n), y)
D <- sqrt( Dx^2 + Dy^2 )


###################################################
### code chunk number 39: a26a
###################################################
options(width=75)


###################################################
### code chunk number 40: a26
###################################################
d <- 1
H <- matrix(1, n, d)
F <- 1
Q <- 1


###################################################
### code chunk number 41: a27 (eval = FALSE)
###################################################
## for( alpha in I(2^(-5:4)) ) { ; #### this will take a while ...
##    Z.tilde <- matrix(NA, tt, n)
##    R <- exp( -alpha*D )
##    for( ii in 1:n ) {
##       cat(ii, " ")
##       bb <- solve( R[ -ii, -ii] ) %*% R[ -ii, ii]
##       xslv <- SS.solve(Z[ , -ii], F=F, Q=Q, H=H[-ii, , drop=FALSE], 
##          R=R[ -ii, -ii], length.out=tt, P0=10^5, beta0=0)
##       Z.hat <- t(  H[-ii, , drop=FALSE] %*% t(xslv$B.apos)  )
##       z.0 <- H[ ii, , drop=FALSE ] %*% t(xslv$B.apos)
##       cov.adj <- (Z[ , -ii] - Z.hat) %*% bb
##       z.tilde <- t(z.0) + cov.adj
##       Z.tilde[ , ii] <- z.tilde[ , 1]
##    } ;
##    rmse <- sqrt( mean( ( Z - Z.tilde )[10:tt, ]^2 ) )
##    cat( alpha, rmse, "\n" )
## }


###################################################
### code chunk number 42: a28
###################################################
alpha <- 1
R <- exp( -alpha*D )

xslv <- SS.solve(Z=Z, F=F, Q=Q, H=H, R=R, length.out=tt, P0=10^5, beta0=0)

x.grid <- seq( min(x)-0.5, max(x)+0.5, length=100 )
y.grid <- seq( min(y)-0.5, max(y)+0.5, length=100 )
xdom <- rep( x.grid, length(y.grid) )
ydom <- rep( y.grid, each=length(x.grid) )

Dx.dom <- tcrossprod(x, rep(1, length(xdom))) - tcrossprod(rep(1, n), xdom)
Dy.dom <- tcrossprod(y, rep(1, length(ydom))) - tcrossprod(rep(1, n), ydom)
D.dom <- t(  sqrt( Dx.dom^2 + Dy.dom^2 )  )

R0 <- exp(-alpha*D.dom)
bb <- solve(R) %*% t(R0)

Hdom <- matrix(1, length(xdom), 1)


###################################################
### code chunk number 43: a29 (eval = FALSE)
###################################################
## for(i in 1:nrow(Z)) {
##    z.hat <- H %*% xslv$B.apos[i, ]
##    z.0 <- Hdom %*% xslv$B.apos[i, ]
##    z.tilde <- z.0 + t( t(Z[i, ] - z.hat) %*% bb )
##    Z.mx <- matrix( z.tilde, length(y.grid), length(x.grid) )
##    image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000))
##    points(x, y, cex=(Z[i ,]-min(Z))/30)
##    map("state", "california", add=TRUE)
##    Sys.sleep(0.1)
## }


###################################################
### code chunk number 44: myFigX
###################################################
png("myFigX.png")
i <- 600
z.hat <- H %*% xslv$B.apos[i, ]
z.0 <- Hdom %*% xslv$B.apos[i, ]
z.tilde <- z.0 + t( t(Z[i, ] - z.hat) %*% bb )
Z.mx <- matrix( z.tilde, length(y.grid), length(x.grid) )
image(x.grid, y.grid, Z.mx, zlim=range(Z), main=xdate[i], col=heat.colors(10000))
points(x, y, cex=(Z[i ,]-min(Z))/30)
map("state", "california", add=TRUE)
dev.off()


