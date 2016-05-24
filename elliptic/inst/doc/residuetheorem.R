### R code from vignette source 'residuetheorem.Rnw'

###################################################
### code chunk number 1: requirepackage
###################################################
require(elliptic,quietly=TRUE)


###################################################
### code chunk number 2: chooseR
###################################################
R <- 400


###################################################
### code chunk number 3: definesemi
###################################################
u1     <- function(x){R*exp(pi*1i*x)}
u1dash <- function(x){R*pi*1i*exp(pi*1i*x)}


###################################################
### code chunk number 4: straightpart
###################################################
u2     <- function(x){R*(2*x-1)}
u2dash <- function(x){R*2}


###################################################
### code chunk number 5: residuetheorem.Rnw:95-96
###################################################
f <- function(z){exp(1i*z)/(1+z^2)}


###################################################
### code chunk number 6: ansapp
###################################################
answer.approximate <-
    integrate.contour(f,u1,u1dash) +
    integrate.contour(f,u2,u2dash) 


###################################################
### code chunk number 7: compareans
###################################################
answer.exact <- pi/exp(1)
abs(answer.approximate - answer.exact)


###################################################
### code chunk number 8: residuetheorem.Rnw:123-124
###################################################
abs(integrate.segments(f,c(-R,R,1i*R))- answer.exact)


###################################################
### code chunk number 9: useabigsquare
###################################################
abs(integrate.segments(f,c(-R,R,R+1i*R, -R+1i*R))- answer.exact)


###################################################
### code chunk number 10: residuetest
###################################################
f <- function(z){sin(z)}
numerical <- residue(f,z0=1,r=1)
exact <- sin(1)
abs(numerical-exact)


