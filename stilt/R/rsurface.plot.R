rsurface.plot <-
function(emul, parind, parvals, tind, n1, n2, zlim=NULL) {



# PRELIMINARIES #!+
m.par     <- dim(emul$Theta.mat)[2] #!+
stopifnot(length(parvals) == m.par, all(parind <= m.par), n1>1, n2>1,
          length(parind) == 2) #!+
if (any(is.finite(parvals[parind]))) cat("WARNING: parvals[parind] should be NA\n")#!+


# CONSTRUCT POINTS WHERE TO PREDICT #!+
# Preliminaries #!+
Theta.mat <- emul$Theta.mat
n.par      <- emul$n
p.par       <- emul$p 
par.min      <- apply(Theta.mat, 2, min)
par.max       <- apply(Theta.mat, 2, max)
par.range      <- par.max - par.min 

# Exceptions #!+
# Where parvals are NA, the comparisons should result in NA
parvals.ok <- (parvals >= par.min) & (parvals <= par.max) 
if (any(!parvals.ok, na.rm=TRUE)) stop("***ERROR*** At least one 'parval' is out of range:)")

# Points where to predict #!+
# x.vec => P1 settings for rows of mu.mat
# y.vec => P2 settings for columns of mu.mat
x.vec     <- seq(par.min[parind[1]], par.max[parind[1]], length.out=n1) 
y.vec     <- seq(par.min[parind[2]], par.max[parind[2]], length.out=n2)


# PREDICT SURFACE #w
mu.mat    <- matrix(NA, nrow=n1, ncol=n2) 
# Mu matrix is filled by columns (y parameter)
for (xind in 1:n1) {
  for (yind in 1:n2) {
     theta.vec          <- parvals 
     theta.vec[parind]  <- c(x.vec[xind], y.vec[yind])
     predict            <- emul.predict(emul, theta.vec) 
     mu.mat[xind, yind] <- predict$mean[tind] 
  }
}


# PLOT SURFACE #!+
# Preliminaries #!+
mar0 <- par("mar")
on.exit(par(mar=mar0)) 
par(mar=c(4,5,5,15)) 

# Specify z-limits #!+
if (is.null(zlim)) zlim=range(mu.mat)

#Plot #!+
myxlab <- paste("Parameter", as.character(parind[1])) 
myylab <- paste("Parameter", as.character(parind[2])) 
filled.contour(x.vec, y.vec, mu.mat, zlim=zlim, color.palette=topo.colors,
               plot.title=title(main="Output", xlab=myxlab, ylab=myylab)) #!+


}
