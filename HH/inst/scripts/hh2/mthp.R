### R code from vignette source '~/WindowsC/HOME/rmh/hh.e2/hh2/mthp.tex'

###################################################
### code chunk number 1: mthp.tex:9-10
###################################################
library(HH)


###################################################
### code chunk number 2: mthp.tex:53-71
###################################################
## hhpdf("line.pdf", height=3, width=6.5)

xyplot(c(-2, 0) ~ c(0, 4), xlim=c(-2.5, 5), ylim=c(-3.5, 1.1),
       xlab="x", ylab=list("y", rot=0),
       scales=list(tck=c(1, .5)),
       pch="+", cex=3, col="#0080ff",
       main=expression(y == -2 + frac(1, 2) ~ x ~ "" == a + b*x),
       par.settings=list(clip=list(panel=FALSE),
                         layout.widths=list(right.padding=7))) +
  layer(panel.abline(a=-2, b=.5, col="#0080ff")) +
  layer(panel.abline(h=0, v=0, col="gray60")) +
  layer(panel.segments(2, -2, c(0,2), c(-2, -1), lty=1)) +
  layer(panel.text(x=c(1, 2.2), y=c(-2.2, -1.5), expression(2, 1))) +
  layer(panel.abline(h=-2, col="gray60", lty=2)) +
  layer(panel.axis("right", at=-2,
        labels=expression(a),
        outside=TRUE, rot=0, line.lty=3, line.col="gray60", line.lwd=2, tck=2))
## hhdev.off()


###################################################
### code chunk number 3: mthp.tex:98-119
###################################################
## hhpdf("parabola.pdf", height=3, width=6.5)
parabola <- data.frame(x=seq(-6, 4, .05))
parabola$y <- with(parabola, (x-3)*(x+5))
xyplot(y ~ x, data=parabola, type="l",
       ylab=list(rot=0),
       scales=list(tck=c(1, .5)),
       main=expression(y == (x+5)*(x-3) ~ ""
                         == 1*x^2 + 2*x - 15 ~ ""
                         == a*x^2 + b*x + c),
       par.settings=list(clip=list(panel=FALSE),
                         layout.widths=list(right.padding=7))) +
  layer(panel.abline(h=0, v=0, col="gray60")) +
  layer(panel.abline(v=c(-5, 3, -1), col="gray60", lty=c(3,3,2))) +
  layer(panel.abline(h=-16, col="gray60", lty=2)) +
  layer(panel.axis("top", at=c(-5, 3, -1),
        labels=expression(root[1], root[2], minimum),
        outside=TRUE, rot=0, line.lty=3, line.col="gray60", line.lwd=2, tck=2)) +
  layer(panel.axis("right", at=-16,
        labels=expression(minimum),
        outside=TRUE, rot=0, line.lty=3, line.col="gray60", line.lwd=2, tck=2))
## hhdev.off()


###################################################
### code chunk number 4: mthp.tex:142-172
###################################################
## hhpdf("ellipse32.pdf", height=4, width=5.5)
a <- 3
b <- 2

theta <- seq(0,(2*pi),length=101)
x <- a * cos(theta)
y <- b * sin(theta)

xyplot( y ~ x,
       aspect=2/3, type="l",
       xlim=c(-3.3,3.3), ylim=c(-2.2,2.2),
       scales=list(at=-3:3, tck=c(1, .5)),
       ylab=list(rot=0),
       main=expression(frac(x^2,3^2) + frac(y^2, 2^2) ==
                       frac(x^2,a^2) + frac(y^2, b^2) ~ "" == 1 ),
       par.settings=list(clip=list(panel=FALSE),
                         layout.widths=list(right.padding=7))) +
         layer(panel.abline(h=0, v=0, col="gray60", lty=3)) +
         layer(panel.abline(h=2, v=3, col="gray60", lty=2)) +
         layer(panel.segments(0, 0, c(0, 3), c(2, 0), lty=1)) +
         layer(panel.axis("top", at=3,
               labels=expression(a),
               outside=TRUE, rot=0, line.lty=3, line.col="gray60",
               line.lwd=2, tck=2)) +
         layer(panel.axis("right", at=2,
              labels=expression(b),
              outside=TRUE, rot=0, line.lty=3, line.col="gray60",
              line.lwd=2, tck=2))

## hhdev.off()


###################################################
### code chunk number 5: mthp.tex:205-219
###################################################
## hhpdf("simultaneous.pdf", height=4, width=5.5)

xyplot( 3 ~ 5, xlim=c(-.2, 9.2), ylim=c(-.8, 4.3),
       xlab="x", ylab=list("y", rot=0), asp=.5,
       scales=list(tck=c(1, .5)),
       pch="+", cex=3, col="#0080ff",
       main=expression(x + y == 8 * "," ~~ 2*x - 3*y == 1)) +
  layer(panel.abline(a=8, b=-1, col="#0080ff")) +
  layer(panel.abline(a=-1/3, b=2/3, col="#0080ff")) +
  layer(panel.abline(h=0, v=0, col="gray60")) +
  layer(panel.text(x=c(7, 3.4, 5.75), y=c(1.5, 1.5, 3), srt=c(-45, 33, 0),
        labels=expression(x+y==8, 2*x-3*y == 1, "(" * 5 * ", " * 3 * ")")))

## hhdev.off()


###################################################
### code chunk number 6: mthp.tex:279-299
###################################################
## hhpdf("explog.pdf", height=4, width=5.5)

tmp <- data.frame(x=seq(-3.1, 6, .05))
tmp$y <- exp(tmp$x)

xyplot(y ~ x, data=tmp, type="l",
       xlim=c(-3.1, 6.1), ylim=c(-3.1, 6.1), aspect=1,
       xlab="x", ylab=list("y", rot=0),
       scales=list(tck=c(1, .5), at=-3:6)) +
  layer(panel.lines(y, x)) +
  layer(panel.abline(h=0, v=0, col="gray60")) +
  layer(panel.abline(h=1, v=1, col="gray60", lty=3)) +
  layer(panel.segments(.5, log(1)-.5*(1/1), 1.5, log(1)+.5*(1/1)))   + ## derivative of ln
##  layer(panel.segments(.5, exp(1)-.5*exp(1), 1.5, exp(1)+.5*exp(1))) + ## derivative of exp at x=1
  layer(panel.segments(-.5, exp(0)-.5*exp(0), .5, exp(0)+.5*exp(0))) + ## derivative of exp at x=0
  layer(panel.text(x=c(2.5, 2.7), y=c(2.5, -1),
                   labels=expression(y == frac(dy, dx) ~ "" == e^x ,
                                     y == ln(x) *~ ";" ~~~ frac(dy, dx) == frac(1, x))))

## hhdev.off()


###################################################
### code chunk number 7: mthp.tex:320-332
###################################################
## hhpdf("asymptote.pdf", height=4, width=5.5)

tmp <- data.frame(x=seq(.05, 10.1, .05))
tmp$y <- 1/tmp$x

xyplot(y ~ x, data=tmp, type="l",
       xlim=c(-.8, 10.1), ylim=c(-.8, 10.1), aspect=1,
       xlab="x", ylab=list("y", rot=0),
       scales=list(tck=c(1, .5), at=0:10)) +
  layer(panel.abline(h=0, v=0, col="gray60"))

## hhdev.off()


###################################################
### code chunk number 8: mthp.tex:461-479
###################################################
## hhpdf("newton.pdf", height=3.5, width=6)

x <- seq(1.75, 2.25, .1)
f <- function(x) x^3 - x - 5
fprime <- function(x) 3*x^2 -1
y <- f(x)

xyplot(y ~ x, type="l", xlim=c(1.75, 2.25)) +
  layer(panel.abline(h=0, v=0, col="gray60")) +
  layer(panel.text(x=c(2+.015, 21/11+.033), y=-.35, label=expression(x[0]==2, x[1]==21/11), col="red")) +
  layer(panel.points(x=c(2, 21/11, 2), y=c(0, 0, 1), pch=c("+","+","x"), col="red", cex=1.5)) +
  layer(panel.segments(2-.15, f(2)-.15*fprime(2), 2+.15, f(2)+.15*fprime(2), col="blue")) +
  layer(panel.segments(2, 0, 2, f(2), col="red")) +
  layer(panel.text(x=2.04, y=.9, labels=expression("(" * 2 * ", " * 1 * ")"))) +
  layer(panel.text(x=c(2.18, 2.05), y=c(2, 2.8),
        labels=expression(y == -21+f*minute(x[0])*x, f(x) == x^3-x-5)))

## hhdev.off()


###################################################
### code chunk number 9: mthp.tex:897-910
###################################################
## hhcapture("qrExample.Rout", '
X <- matrix(c(1,3,6,4,2,3,8,6,4,5,3,2), 4, 3)
X
crossprod(X)

## use the efficient calculation
X.qr <- qr(X)
qr.Q(X.qr) ## display q
qr.R(X.qr) ## display r
zapsmall(crossprod(qr.Q(X.qr)))  ## identity
crossprod(qr.R(X.qr))    ## reproduce crossprod(X)
qr.X(X.qr) ## reproduce X
## ')


###################################################
### code chunk number 10: mthp.tex:991-1023
###################################################
hhcode("mgs.R", '
## modified Gram-Schmidt orthogonalization

mgs <- function(x) {
  ## modified Gram-Schmidt orthogonalization

  ## this is an expository algorithm
  ## this is not an efficient computing algorithm

  ## q[,j] is the normalized residual from the least squares fit of
  ## x[,j] on the preceding normalized columns q[,1:(j-1)]

  n <- nrow(x)
  m <- ncol(x)

  q <- x
  r <- matrix(0, m, m)

  for (i in 1:m) {
    r[i,i] <- sqrt(sum(q[,i]^2)) ## length of q[,i]
    q[,i] <- q[,i] / r[i,i]      ## normalize q[,i]

    if (i < m) {   ## if we still have columns to go
      for (j in (i+1):m) {
        r[i,j] <- sum(q[,i] * q[,j]) ## length of projection of q[,j] on q[,i]
        q[,j] <- q[,j] - q[,i] * r[i,j] ## remove projection of q[,j] on q[,i]
      }
    }
  }
  list(q=q, r=r)
}
## ')


###################################################
### code chunk number 11: mthp.tex:1036-1053
###################################################
hhcode("mgs-example.R", '
X <- matrix(c(1,3,6,4,2,3,8,6,4,5,3,2), 4, 3)
X
crossprod(X)

## use the expository function defined in the previous chunk

X.mgs <- mgs(X)
X.mgs  ## q is orthogonal, r is upper triangular
## These are identical to the results of qr(X)
## up to the sign of the columns of q and the rows of r.

zapsmall(crossprod(X.mgs$q))  ## identity
crossprod(X.mgs$r)    ## reproduces crossprod(X)

X.mgs$q %*% X.mgs$r   ## reproduces X
## ')


###################################################
### code chunk number 12: mthp.tex:1082-1089
###################################################
## hhcapture("cholExample.Rout", '
X <- matrix(c(1,3,6,4,2,3,8,6,4,5,3,2), 4, 3)
M <- crossprod(X)
M
chol(M)
crossprod(chol(M)) ## reproduce M
## ')


###################################################
### code chunk number 13: mthp.tex:1133-1144
###################################################
## hhcapture("projExample.Rout", '
X <- matrix(c(3,1,0, 1,2,0, 0,0,0), 3, 3)
P <- cbind(qr.Q(qr(X))[, 1:2], 0)
P
crossprod(P)
y <- matrix(1:3)
y
P %*% y
sqrt(sum(y[1:2,]^2))
sqrt(sum((P %*% y)^2))
## ')


###################################################
### code chunk number 14: mthp.tex:1187-1220
###################################################
## hhpdf("rotate.pdf", height=4, width=4)

u <- matrix(c(1.5, 1))
theta <- -pi/7 ## radians
M <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), 2, 2)
M
v <- M %*% u
v
axes <- matrix(c(1, 0, 0, 1), 2, 2)
M.axes <- M %*% axes

xyplot(c(0, 1) ~ c(1, 0), xlim=c(-2.25, 2.25), ylim=c(-2.25, 2.25), aspect=1,
       pch=19, col="red", xlab="x", ylab="y") +
  layer(panel.abline(h=0, v=0, col="gray60")) +
  layer(panel.segments(axes[1,], axes[2,], 0, 0, col="red")) +
  layer(panel.points(u[1], u[2], pch=19, col="blue")) +
  layer(panel.text(u[1]+.2, u[2], label="u", col="blue")) +
  layer(panel.segments(u[1], u[2], 0, 0, col="blue")) +
  layer(panel.points(v[1], v[2], pch=19, col="orange")) +
  layer(panel.text(v[1]+.15, v[2], label="v", col="orange")) +
  layer(panel.segments(v[1], v[2], 0, 0, col="orange")) +
  layer(panel.segments(M.axes[1,], M.axes[2,], 0, 0, col="green", lwd=2)) +
  layer(panel.points(M.axes[1,], M.axes[2,], pch=19, col="green", lwd=2)) +
  layer(panel.abline(a=0, b=tan(-theta),      col="green", lty=2)) +
  layer(panel.abline(a=0, b=tan(-theta-pi/2), col="green", lty=2)) +
  layer(panel.text(.40, .07, labels=expression(theta)))

trellis.focus(highlight=FALSE)
grid.curve(.3, 0, .3*cos(theta), .3*sin(-theta), angle=(theta/pi)*180,
           gp=gpar(col="purple"), default.units="native")
trellis.unfocus()

## hhdev.off()


###################################################
### code chunk number 15: mthp.tex:1302-1313
###################################################
## hhcapture("eigen2111.Rout", '
V <- matrix(c(2, 1, 1, 1), 2, 2)
V
eV <- eigen(V)
eV
sqrt(eV$val)  ## semimajor and semiminor axis lengths
## angle of axes in radians
atan(c(eV$vec[2,1]/eV$vec[1,1], eV$vec[2,2]/eV$vec[1,2]))
## = -pi/2 ## right angle
diff(atan(c(eV$vec[2,1]/eV$vec[1,1], eV$vec[2,2]/eV$vec[1,2])))
## ')


###################################################
### code chunk number 16: mthp.tex:1326-1351
###################################################
## hhpdf("ellipse2111.pdf", width=5.5, height=4)
theta <- seq(0,(2*pi),length=101)
x <- sqrt(eV$val[1]) * cos(theta)
y <- sqrt(eV$val[2]) * sin(theta)

newxy <- cbind(x,y) %*% t(eV$vec)

xyplot( newxy[,2] ~ newxy[,1], type="l",
       xlab="x", ylab=list("y", rot=0),
       scales=list(tck=c(1, .5), at=seq(-1.5,1.5,.5)),
       aspect=2/3, col="mediumblue",
       xlim=c(-1.65, 1.65), ylim=c(-1.1, 1.1)) +
         layer(panel.abline(h=0, v=0, col="gray60")) +
         layer(panel.abline(a=0, b=eV$vec[2,1]/eV$vec[1,1], lty=3, col="gray30")) +
         layer(panel.abline(a=0, b=eV$vec[2,2]/eV$vec[1,2], lty=3, col="gray30")) +
         layer(panel.segments(0, 0, eV$vec[1,1:2]*sqrt(eV$val), eV$vec[2,1:2]*sqrt(eV$val), lty=2, lwd=2)) +
         layer(panel.segments(sum(eV$vec[1,1:2])/15,
                              sum(eV$vec[2,1:2])/15,
                              eV$vec[1,1:2]/15,
                              eV$vec[2,1:2]/15, lty=1)) +
         layer(panel.text(x=c(-.7, .3), y=c(-.3, -.25), labels=expression(sqrt(lambda[1]), sqrt(lambda[2])))) +
         layer(panel.text(x= eV$vec[1,] + .06, y= eV$vec[2,] - .05, labels=expression(xi[1], xi[2]))) +
xyplot(sin(theta) ~ cos(theta), type="l", lty=3, lwd=2, col="gray30") ## unit circle

## hhdev.off()


###################################################
### code chunk number 17: mthp.tex:1460-1469
###################################################
## hhcapture("svd1.Rout", '
M <- matrix(c(1,3,6,4,2,3,8,6,4,5,3,2), 4, 3)
M
M.svd <- svd(M)
M.svd
zapsmall(crossprod(M.svd$u))
zapsmall(crossprod(M.svd$v))
M.svd$u %*% diag(M.svd$d) %*% t(M.svd$v)
## ')


###################################################
### code chunk number 18: mthp.tex:1493-1498
###################################################
## hhcapture("svd2.Rout", '
eigen(tcrossprod(M))
eigen(crossprod(M))
M.svd$d^2
## ')


###################################################
### code chunk number 19: mthp.tex:1527-1537
###################################################
## hhcapture("ginv.Rout", '
M <- matrix(c(1,3,6,4,2,3,8,6,4,5,3,2), 4, 3)
M
library(MASS)
Mi <- ginv(M)
Mi
zapsmall(eigen(M %*% Mi)$value)
zapsmall(Mi %*% M)
M.svd$v %*% diag(1/M.svd$d) %*% t(M.svd$u)
## ')


