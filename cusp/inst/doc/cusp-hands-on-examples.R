### R code from vignette source 'cusp-hands-on-examples.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: cusp-hands-on-examples.Rnw:41-42
###################################################
library(cusp)


###################################################
### code chunk number 2: figZeemanMachineSchema
###################################################
plot.new()
plot.window(c(-1.05,1.05),c(-1.1,1))
rect(-.25,-1,.25,-.5, col=rgb(.3,.6,.5,.1),border=NA)
text(.25,-.75, "control plane", pos=4, offset=1)
arrows(-.25,-1,-.25,-.75, len=.1, lwd=.5)            # y arrow
text(-.25,-.875,"y",pos=2)
arrows(-.25,-1,   0,  -1, len=.1, lwd=.5)            # x arrow
text(-.125,-1,"x",pos=1)
points(0, 0.75, pch=20, cex=2)                       # fixed point
text(0,.75,"fixed point",pos=4, offset=1)
abline(v=0)                                          # central axis
text(0,1,"central axis",pos=4, offset=0.2)
y = seq(0,2*pi,len=300)
lines(.3*cos(y),.3*sin(y))                           # circle
arrows(.3*cos(y[30]), .3*sin(y[30]),  .0,  0.75, len=0, lwd=3, col="gray") #coord1
arrows(.3*cos(y[30]), .3*sin(y[30]), .15, -0.75, len=0, lwd=3, col="gray") #coord2
points(.15, -.75, pch=20, col=2)                      # control point
points(0,0,pch=20) # circle center
points(.3*cos(y[30]),.3*sin(y[30]), pch=1)           # strap point
text(.3*cos(y[30]),.3*sin(y[30]), "strap point", pos=4, offset=.5)
arrows(.3*cos(y[30]), .3*sin(y[30]), 0, .3*sin(y[30]), col=4, len=.1)
text(0.5*.3*cos(y[30]),.3*sin(y[30]), "z", pos=1, offset=.5, col=4)


###################################################
### code chunk number 3: figZeemanMachineSchema
###################################################
plot.new()
plot.window(c(-1.05,1.05),c(-1.1,1))
rect(-.25,-1,.25,-.5, col=rgb(.3,.6,.5,.1),border=NA)
text(.25,-.75, "control plane", pos=4, offset=1)
arrows(-.25,-1,-.25,-.75, len=.1, lwd=.5)            # y arrow
text(-.25,-.875,"y",pos=2)
arrows(-.25,-1,   0,  -1, len=.1, lwd=.5)            # x arrow
text(-.125,-1,"x",pos=1)
points(0, 0.75, pch=20, cex=2)                       # fixed point
text(0,.75,"fixed point",pos=4, offset=1)
abline(v=0)                                          # central axis
text(0,1,"central axis",pos=4, offset=0.2)
y = seq(0,2*pi,len=300)
lines(.3*cos(y),.3*sin(y))                           # circle
arrows(.3*cos(y[30]), .3*sin(y[30]),  .0,  0.75, len=0, lwd=3, col="gray") #coord1
arrows(.3*cos(y[30]), .3*sin(y[30]), .15, -0.75, len=0, lwd=3, col="gray") #coord2
points(.15, -.75, pch=20, col=2)                      # control point
points(0,0,pch=20) # circle center
points(.3*cos(y[30]),.3*sin(y[30]), pch=1)           # strap point
text(.3*cos(y[30]),.3*sin(y[30]), "strap point", pos=4, offset=.5)
arrows(.3*cos(y[30]), .3*sin(y[30]), 0, .3*sin(y[30]), col=4, len=.1)
text(0.5*.3*cos(y[30]),.3*sin(y[30]), "z", pos=1, offset=.5, col=4)


###################################################
### code chunk number 4: zeeman1
###################################################
data(zeeman1)
nrow(zeeman1)
head(zeeman1)


###################################################
### code chunk number 5: fit1zeeman1
###################################################
fit1.1 = cusp(y~z, alpha~x+y, beta~x+y, zeeman1)
summary(fit1.1)


###################################################
### code chunk number 6: fit2zeeman1
###################################################
fit1.2 = cusp(y~z-1, alpha~x-1, beta~y, zeeman1)
summary(fit1.2) # compare with logistic fit as well


###################################################
### code chunk number 7: fit3zeeman1
###################################################
fit1.3 = cusp(y~z-1, alpha~x, beta~y, zeeman1)
(sf1.3 <- summary(fit1.3, logist=TRUE))


###################################################
### code chunk number 8: figfit3zeeman1
###################################################
plot(fit1.3)


###################################################
### code chunk number 9: figfit3zeeman1
###################################################
plot(fit1.3)


###################################################
### code chunk number 10: cusp-hands-on-examples.Rnw:164-167
###################################################
data(zeeman2)
fit2.1 = cusp(y~z, alpha~x+y, beta~x+y, zeeman2)
summary(fit2.1)


###################################################
### code chunk number 11: cusp-hands-on-examples.Rnw:170-172
###################################################
fit2.2 = cusp(y~z-1, alpha~x-1, beta~y, zeeman2)
summary(fit2.2)


###################################################
### code chunk number 12: cusp-hands-on-examples.Rnw:175-177
###################################################
fit2.3 = cusp(y~z-1, alpha~x, beta~y, zeeman2)
fit2.3


###################################################
### code chunk number 13: cusp-hands-on-examples.Rnw:183-186
###################################################
data(zeeman3)
fit3.1 <- cusp(y~z, alpha~x+y, beta~x+y, zeeman3)
summary(fit3.1)


###################################################
### code chunk number 14: cusp-hands-on-examples.Rnw:189-191
###################################################
fit3.2 <- cusp(y~z-1, alpha~x-1, beta~y, zeeman3)
summary(fit3.2)


###################################################
### code chunk number 15: hist
###################################################
set.seed(423)
alpha = 0.25
beta = 2
n = 1000
y = rcusp(n, alpha, beta)
hist(y,80,freq=FALSE)
curve(dcusp(x, 0.25, 2), min(y)-1, max(y)+1, add=TRUE, col=2)


###################################################
### code chunk number 16: fighist
###################################################
set.seed(423)
alpha = 0.25
beta = 2
n = 1000
y = rcusp(n, alpha, beta)
hist(y,80,freq=FALSE)
curve(dcusp(x, 0.25, 2), min(y)-1, max(y)+1, add=TRUE, col=2)


###################################################
### code chunk number 17: cusp-hands-on-examples.Rnw:227-228
###################################################
cusp(y~y-1, alpha~1, beta~1)


###################################################
### code chunk number 18: cusp-hands-on-examples.Rnw:232-241
###################################################
set.seed(423)
x1 = runif(150)
x2 = runif(150)
a = c(-2, 4)
b = c(-1, 4)
alpha = a[1] + a[2]*x1
beta  = b[1] + b[2]*x2
z = Vectorize(rcusp)(1, alpha, beta)
data <- data.frame(x1, x2, z)


###################################################
### code chunk number 19: cusp-hands-on-examples.Rnw:244-245
###################################################
fit <- cusp(y ~ z, alpha ~ x1+x2, beta ~ x1+x2, data)


###################################################
### code chunk number 20: generate-measurement-error
###################################################
set.seed(423)
g = expand.grid(seq(-3,3,len=15), seq(-3,3,len=15))
a = g[,1]; b = g[,2]; idx=cbind(sample(c(1,3),length(a),TRUE), seq(along=a));
s = Vectorize(cusp.extrema)(a,b)[idx]
y = s + rnorm(length(s),,.3)


###################################################
### code chunk number 21: 3Dscatter
###################################################
if(require(plot3D)){
	scatter3D(a, b, y, theta=200, phi=10, zlim=c(-4,4));
  # you can try require(rgl); rgl.points(a, b, y) instead
}


###################################################
### code chunk number 22: fig1
###################################################
if(require(plot3D)){
	scatter3D(a, b, y, theta=200, phi=10, zlim=c(-4,4));
  # you can try require(rgl); rgl.points(a, b, y) instead
}


###################################################
### code chunk number 23: cusp-hands-on-examples.Rnw:279-281
###################################################
fit <- cusp(y~y, alpha~a, beta~b)
summary(fit)


