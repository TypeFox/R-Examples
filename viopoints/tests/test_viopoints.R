library(viopoints)

methods(viopoints)

viopoints(iris[1])
viopoints(iris[1:4])

# groups
at <- viopoints(iris[1], groups=iris[5], col=2:4)
at
at <- viopoints(iris[1], groups=iris[5], at=1, col=2:4)
at

# data.frame
viopoints(iris[1:4], groups=iris[5], at=1:4, col=2:4)
viopoints(iris[,1:4], groups=iris[5], at=1:4, col=2:4)
viopoints(iris[1:4], groups=iris[,5], at=1:4, col=2:4)
viopoints(iris[,1:4], groups=iris[,5], at=1:4, col=2:4)

# list
viopoints(list(unif=runif(300,-1,1), norm=rnorm(500), exp=rexp(100)), col=2:4)

# vector, matrix, factor
viopoints(as.vector(iris[,1]), groups=iris[,5], at=1:4, col=2:4)
viopoints(as.vector(iris[,1]), groups=as.factor(iris[,5]), at=1:4, col=2:4)
viopoints(as.matrix(iris[,1:4]), groups=iris[,5], at=1:4, col=2:4)
viopoints(as.matrix(iris[,1:4]), groups=as.factor(iris[,5]), at=1:4, col=2:4)
viopoints(iris[1:4], groups=as.factor(iris[,5]), at=1:4, col=2:4)
viopoints(iris[,1:4], groups=as.factor(iris[,5]), at=1:4, col=2:4)

# formula
viopoints(Sepal.Length ~ Species, data=iris, col=2:4)
viopoints(Sepal.Length ~ 1, data=iris, col=2:4)
viopoints(as.matrix(iris[,1:4]) ~ iris[,5], at=1:4, col=2:4)
viopoints(as.matrix(iris[,1:4]) ~ 1)
viopoints(Sepal.Width ~ Species + cut(Sepal.Width, 2), data=iris, col=rainbow(6, v=0.8))
viopoints(Sepal.Width ~ Species + cut(Sepal.Width, 2), data=iris, at=1:3, col=rainbow(6, v=0.8))

x <- stats::rnorm(100)

# method
method <- c("overplot", "jitter", "stack", "violin")
viopoints(x, method="overplot", xlim=c(0,4), at=0.5, glab="method", axis.par=list(at=seq(0.5,3.5), labels=method))
viopoints(x, method="jitter", at=1.5, add=TRUE)
viopoints(round(x,1), method="stack", at=2.5, add=TRUE)
viopoints(x, method="violin", at=3.5, add=TRUE)

viopoints(x, method="overplot", horizontal=TRUE, ylim=c(0,4), at=0.5, glab="method", axis.par=list(at=seq(0.5,3.5), labels=method))
viopoints(x, method="jitter", horizontal=TRUE, at=1.5, add=TRUE)
viopoints(round(x,1), method="stack", horizontal=TRUE, at=2.5, add=TRUE)
viopoints(x, method="violin", horizontal=TRUE, at=3.5, add=TRUE)

# amount, offset
viopoints(round(x,1), method="stack", side="p", at=0.5, xlim=c(0,4), amount=0.2, glab="amount", main="stack", axis.par=list(at=seq(0.5,3.5), labels=c("0.3","0.5","NA","NA,0.5")))
viopoints(round(x,1), method="stack", side="p", at=1.5, amount=0.3, add=TRUE) 
viopoints(round(x,1), method="stack", side="p", at=2.5, col=3, amount=NA, add=TRUE) 
viopoints(round(x,1), method="stack", side="p", at=3.5, col=3, amount=NA, offset=0.5, add=TRUE) 

viopoints(round(x,1), horizontal=TRUE, method="stack", side="p", at=0.5, ylim=c(0,4), amount=0.2, glab="amount", main="stack", axis.par=list(at=seq(0.5,3.5), labels=c("0.3","0.5","NA","NA,0.5")))
viopoints(round(x,1), horizontal=TRUE, method="stack", side="p", at=1.5, amount=0.3, add=TRUE) 
viopoints(round(x,1), horizontal=TRUE, method="stack", side="p", at=2.5, col=3, amount=NA, add=TRUE) 
viopoints(round(x,1), horizontal=TRUE, method="stack", side="p", at=3.5, col=3, amount=NA, offset=0.5, add=TRUE) 

# side
side <- c("negative", "both", "positive")
viopoints(x, x, x, side=side, main="rnorm(100)", glab="side", group.names=side)
viopoints(x, x, x, horizontal=TRUE, side=side, main="rnorm(100)", glab="side", group.names=side)

# main, sub, xlab, ylab, dlab, glab
viopoints(iris[1:4], main="main", sub="sub", xlab="xlab", ylab="ylab")
viopoints(iris[1:4], horizontal=TRUE, main="main", sub="sub", xlab="xlab", ylab="ylab")
viopoints(iris[1:4], main="main", sub="sub", dlab="dlab", glab="glab")
viopoints(iris[1:4], horizontal=TRUE, main="main", sub="sub", dlab="dlab", glab="glab")

# lines, line.lty, line.lwd, line.col 
viopoints(iris[,1:4], groups=iris[,5], col=2:4, lines=TRUE)
viopoints(iris[,1:4], groups=iris[,5], col=2:4, lines=TRUE, line.col=c("pink", "lightgreen", "skyblue"))
viopoints(iris[,1:4], groups=iris[,5], at=1:4, col=2:4, lines=TRUE, line.col=c("pink", "lightgreen", "skyblue"))
viopoints(iris[,1:4], groups=iris[,5], at=1:4, col=2:4, lines=TRUE, line.lty=3, line.lwd=2, line.col=c("pink", "lightgreen", "skyblue"))
viopoints(iris[,1:4], groups=iris[,5], at=1:4, horizontal=TRUE, col=2:4, lines=TRUE, line.col=c("pink", "lightgreen", "skyblue"))

# points
viopoints(iris[,1:4], groups=iris[,5], at=1:4, points=F, lines=TRUE, line.col=c("pink", "lightgreen", "skyblue"))

# pch, cex, col, bg
viopoints(iris[,1:4], groups=iris[,5], at=1:4, pch=21:24, cex=1.5, col=2:4, bg=c("pink", "lightgreen", "skyblue"))

# add
boxplot(Sepal.Length ~ Species, data=iris)
viopoints(Sepal.Length ~ Species, data=iris, col=2:4, add=TRUE)

boxplot(iris[1:4])
viopoints(iris[,1:4], groups=iris[,5], at=1:4, col=2:4, lines=TRUE, line.col=c("pink", "lightgreen", "skyblue"), add=TRUE)

# log, ylim
viopoints(decrease ~ treatment, data = OrchardSprays, col=rainbow(8, v=0.8), log = "y", ylim = c(1,200))

# group.names, xlim, jitter
viopoints(list(runif(300,-1,1), rnorm(500), rexp(100)), xlim=c(-5,5), group.names=c("unif", "norm", "exp"), horizontal=T, jitter=0.3, col=2:4)

# axis.par, jitter
viopoints(as.matrix(iris[1:4]) ~ Species, data=iris, jitter=0.3, col=2:4, axis.par=list(at=1:12, labels=paste(rep(substr(levels(factor(iris[,5])),1,3), each=4), c("SL","SW","PL","PW"), sep="."), cex.axis=0.8, las=3))

# density.par
viopoints(iris[3], density.par=list(na.rm=TRUE, bw="nrd0"), main="nrd0")
viopoints(iris[3], density.par=list(na.rm=TRUE, bw="SJ-ste"), main="SJ-ste")

# frame.plot, axes
viopoints(iris[1:4], axes=FALSE)
viopoints(iris[1:4], frame.plot=FALSE, axes=TRUE)
viopoints(iris[1:4], frame.plot=TRUE, axes=FALSE)

# NA
na_iris <- iris
na_iris[1,1] <- NA
na_iris[2,5] <- NA
na_iris[3,1] <- NA
na_iris[3,5] <- NA
viopoints(na_iris[1])
viopoints(na_iris[1], groups=na_iris[5])
viopoints(Sepal.Length ~ 1, data=na_iris)
viopoints(Sepal.Length ~ Species, data=na_iris)
viopoints(na_iris[1], groups=na_iris[5], na.group=TRUE)
viopoints(Sepal.Length ~ Species, data=na_iris, na.action=na.pass, na.group=TRUE)
viopoints(na_iris[2], groups=na_iris[5], na.group=TRUE)
viopoints(Sepal.Width ~ Species, data=na_iris, na.action=na.pass, na.group=TRUE)

for (method in c("overplot", "jitter", "stack", "violin")) {
  viopoints(na_iris[1:4], method=method, main=method)
  viopoints(na_iris[3], groups=na_iris[5], method=method, col=2:4, main=method)
  viopoints(na_iris[3], groups=na_iris[5], na.group=TRUE, method=method, col=2:4, main=method)
}


