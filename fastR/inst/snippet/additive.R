x1 = rep(rep(1:4, each=4),4);  x2 = rep(rep(1:4, times=4),4)
fit = 1 + 2 * x1 - 4 * x2; e = rnorm(4*4*4,sd=1)
plot1 <- xyplot(fit~x2,groups=x1,type='r',
      main=expression(paste("Model fits for fixed values of ",x[1])))
plot2 <- xyplot(fit~x1,groups=x2,type='r',
      main=expression(paste("Model fits for fixed values of ",x[2])))
plot3 <- xyplot(I(fit+e)~x2,groups=x1,type=c('p','r'),pch=16,
      main=expression(paste("Simulated data grouped by ",x[1])))
plot4 <- xyplot(I(fit+e)~x1,groups=x2,type=c('p','r'),pch=16,
      main=expression(paste("Simulated data grouped by ",x[2])))
plot5 <- xyplot(I(fit+4*e)~x2,groups=x1,type=c('p','r'),pch=16,
      main=expression(paste("Simulated data grouped by ",x[1], 
                "(larger ", sigma,")")))
plot6 <- xyplot(I(fit+4*e)~x1,groups=x2,type=c('p','r'),pch=16,
      main=expression(paste("Simulated data grouped by ",x[2], 
                "(larger ", sigma,")")))
