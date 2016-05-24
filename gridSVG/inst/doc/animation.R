### R code from vignette source 'animation.Rnw'

###################################################
### code chunk number 1: animation.Rnw:28-30
###################################################
library(grid)
library(gridSVG)


###################################################
### code chunk number 2: circle
###################################################
grid.circle(.1, .5, r=.1, gp=gpar(fill="black"),
            name="circle")
grid.animate("circle", x=c(.1, .9))
grid.export("animCircle.svg")


###################################################
### code chunk number 3: polylinesrc (eval = FALSE)
###################################################
## grid.rect()
## grid.polyline(c(.2, .3, .4, .6, .7, .8),
##               c(.7, .5, .7, .3, .5, .3),
##               id=rep(1:2, each=3),
##               gp=gpar(lwd=5),
##               name="polyline")


###################################################
### code chunk number 4: polyline
###################################################
grid.rect()
grid.polyline(c(.2, .3, .4, .6, .7, .8),
              c(.7, .5, .7, .3, .5, .3),
              id=rep(1:2, each=3),
              gp=gpar(lwd=5),
              name="polyline")


###################################################
### code chunk number 5: animvalue
###################################################
polylineY <- animUnit(unit(c(.7, .5, .7, .3, .5, .3,
                             .3, .5, .3, .7, .5, .7,
                             .7, .5, .7, .3, .5, .3),
                           unit="npc"),
                      timeid=rep(1:3, each=6),
                      id=rep(rep(1:2, each=3), 3))


###################################################
### code chunk number 6: animation.Rnw:83-84
###################################################
polylineY


###################################################
### code chunk number 7: animation.Rnw:116-117
###################################################
animUnit(unit(1:4, "cm"))


###################################################
### code chunk number 8: animation.Rnw:140-141
###################################################
animUnit(unit(1:4, "cm"), id=rep(1:2, 2))


###################################################
### code chunk number 9: animation.Rnw:152-153
###################################################
animUnit(unit(1:12, "cm"), timeid=rep(1:2, 6))


###################################################
### code chunk number 10: animation.Rnw:163-165
###################################################
animUnit(unit(1:12, "cm"), 
         timeid=rep(1:2, 6), id=rep(1:2, each=6))


###################################################
### code chunk number 11: animpolylinesrc (eval = FALSE)
###################################################
## grid.animate("polyline", y=polylineY, rep=TRUE)


###################################################
### code chunk number 12: animpolyline
###################################################
grid.newpage()
grid.rect()
grid.polyline(c(.2, .3, .4, .6, .7, .8),
              c(.7, .5, .7, .3, .5, .3),
              id=rep(1:2, each=3),
              gp=gpar(lwd=5),
              name="polyline")
polylineY <- animUnit(unit(c(.7, .5, .7, .3, .5, .3,
                             .3, .5, .3, .7, .5, .7,
                             .7, .5, .7, .3, .5, .3),
                           unit="npc"),
                      timeid=rep(1:3, each=6),
                      id=rep(rep(1:2, each=3), 3))
grid.animate("polyline", y=polylineY, rep=TRUE)
grid.export("animPolyline.svg")


###################################################
### code chunk number 13: animation.Rnw:197-198
###################################################
as.animUnit(1:4, unit="cm")


###################################################
### code chunk number 14: animation.Rnw:208-209 (eval = FALSE)
###################################################
## grid.animate("circle", x=c(.1, .9))


###################################################
### code chunk number 15: animation.Rnw:217-219
###################################################
m <- matrix(1:6, ncol=2)
m


###################################################
### code chunk number 16: animation.Rnw:222-223
###################################################
as.animUnit(m, unit="cm")


###################################################
### code chunk number 17: animation.Rnw:231-232
###################################################
as.animUnit(m, unit="cm", multVal=TRUE)


###################################################
### code chunk number 18: animation.Rnw:246-248
###################################################
l <- list(unit(1:3, "cm"), unit(4:6, "cm"))
l


###################################################
### code chunk number 19: animation.Rnw:251-252
###################################################
as.animUnit(l)


###################################################
### code chunk number 20: animation.Rnw:260-261
###################################################
as.animUnit(l, multVal=TRUE)


