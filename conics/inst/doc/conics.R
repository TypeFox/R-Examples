### R code from vignette source 'conics.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: conics.Rnw:341-342
###################################################
library(conics)


###################################################
### code chunk number 2: conics.Rnw:362-363
###################################################
v <- c(2,2,2,-20,-28,10)


###################################################
### code chunk number 3: conics.Rnw:366-367
###################################################
A <- conicMatrix(v)


###################################################
### code chunk number 4: conics.Rnw:373-375
###################################################
conicCenter(v)
conicAxes(v)


###################################################
### code chunk number 5: conics.Rnw:378-380 (eval = FALSE)
###################################################
## conicCenter(A)
## conicAxes(A)


###################################################
### code chunk number 6: conics.Rnw:386-387
###################################################
conicPlot(v, main="conicPlot(v)", xlab="", ylab="")


###################################################
### code chunk number 7: conics.Rnw:480-481
###################################################
conicPlot(v, center=T, sym.axes=T)


###################################################
### code chunk number 8: conics.Rnw:490-496
###################################################
v <- c(2,2,2,-20,-28,10)
conicPlot(v, center=T, lty=1)
v[6] <- 30
conicPlot(v, add=T, lty=2)
v[6] <- 50
conicPlot(v, add=T, lty=4)


###################################################
### code chunk number 9: conics.Rnw:505-507
###################################################
v <- c(2,2,-2,-20,20,10)
conicPlot(v, xlim=c(-10,10), ylim=c(-5,18))


###################################################
### code chunk number 10: conics.Rnw:516-519
###################################################
conicPlot(v, asymptotes=T, sym.axes=T,
as.col="red", as.lty=2, ax.col="blue", ax.lty=4, 
xlim=c(-10,10), ylim=c(-5,18))


###################################################
### code chunk number 11: conics.Rnw:528-531
###################################################
conicPlot(v, asymptotes=T, sym.axes=T, 
xlim=c(-10,10), ylim=c(-5,18),
asp=1, col="blue", main="Hyperbola", bty="n")


###################################################
### code chunk number 12: conics.Rnw:583-585
###################################################
v <- c(2,2,-2,-20,20,10)
res <- conicPlot(v)


###################################################
### code chunk number 13: conics.Rnw:589-595
###################################################
v <- c(-4,0,1,0,0,1)
cp <- conicPlot(v, sym.axes=TRUE, asymptote=TRUE, asp=1, 
                   ax.lty=2, as.col="gray")
points(cp$foci$x,cp$foci$y,col="red",pch=19)
text(cp$foci$x,cp$foci$y+0.1,paste("F",2:1,sep=""))
points(cp$vertices$x,cp$vertices$y,col="blue",pch=19)


