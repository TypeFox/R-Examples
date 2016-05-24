# phaseportdat produces a phase portrait based on the best models (differential equations)
# with data trajectories of specified entities (e.g. countries)
# for x and y, the phase portrait is the visual simultaneous solution of the two differential equations 

# to call: phaseportdat(datap, datap$logGDP, datap$EmanzV, seq(0, 1, by = 0.1), seq(0, 1, by = 0.1), f <- function(t,Y=c()) rbind(0.0012/Y[1]^2, + 0.0071*Y[1]^3), 1, 2, 4, 5, 7, 9)
# with dx  = + 0.0012 /x^2 and dy = + 0.0071 x^3 

phaseportdat <- function(dataset, xv, yv, rangeX, rangeY, f, entidx1, entidx2, 
                      entidx3, entidx4, entidx5, entidx6) 
{  
  procdata <- preprocess_data(2, xv, yv)
  xwide <- procdata$xwide
  ywide <- procdata$ywide
#   
#   # boundaries of the system, here under assumption that both variables are scaled 0-1
#   # the visualisation is optimized for this scaling, in case the user would like to use
#   # the original scale, he/she needs to adjust the code
  tmpx <- rangeX
  tmpy <- rangeY
  xmin=min(tmpx)
  xmax=max(tmpx)
  ymin=min(tmpy)
  ymax=max(tmpy)
  rgrid = 25

  y1 <- linspace(xmin, xmax, rgrid)
  y2 <- linspace(ymin, ymax, rgrid)
  
  # scaling arrows of the vectorfield
  qscale <- range(y1)/range(y2)
  mgrid <- meshgrid(y1, y2)
  x <- mgrid[[1]] 
  y <- mgrid[[2]] 
  
  u <- matrix(0, rgrid, rgrid) 
  v <- matrix(0, rgrid, rgrid)
  
  # compute derivates at each point (every possible intial conditions)
  t = 0
  for (i in 1:length(x))
  {
    Yprime <- f(t, rbind(x[i], y[i]))
    u[i] <- Yprime[1]
    v[i] <- Yprime[2]
  }

  # phase portrait with highlighted data trajectories
  dev.set(1)
  postscript("ppdat.eps", horizontal=FALSE, width=5, height=5,
           onefile=FALSE, paper="special", family="ComputerModern")
  
  # setting a plot 
  plot(x, y, col="white", xlab = "X-Variable", ylab = "Y-Variable")
  grid(col="white")
  # velocity plot
  quiver(x, y, u, v, scale=0.5, angle = 10, length = 0.05, lwd = 1.4, col = 'gray') 
  
  # highligh selected entities' tranjectories  
  matplot(xwide[entidx1,], ywide[entidx1,], type='l', ldw=2, col = 'blue', add=TRUE)
  matplot(xwide[entidx2,], ywide[entidx2,], type='l', ldw=2, col = 'darkgreen', add=TRUE)
  matplot(xwide[entidx3,], ywide[entidx3,], type='l', ldw=2, col = 'red', add=TRUE)
  matplot(xwide[entidx4,], ywide[entidx4,], type='l', ldw=2, col = 'cyan', add=TRUE)
  matplot(xwide[entidx5,], ywide[entidx5,], type='l', ldw=2, col = 'magenta', add=TRUE)
  matplot(xwide[entidx6,], ywide[entidx6,], type='l', ldw=2, col = 'black', add=TRUE)
  
  # adding markers (circles) for marking initial conditions (default)
  points(xwide[entidx1,1], ywide[entidx1,1], pch = 20, cex=1.2, col = 'blue')
  points(xwide[entidx2,1], ywide[entidx2,1], pch = 20, cex=1.2, col = 'darkgreen')
  points(xwide[entidx3,1], ywide[entidx3,1], pch = 20, cex=1.2, col = 'red')
  points(xwide[entidx4,1], ywide[entidx4,1], pch = 20, cex=1.2, col = 'cyan')
  points(xwide[entidx5,1], ywide[entidx5,1], pch = 20, cex=1.2, col = 'magenta')
  points(xwide[entidx6,1], ywide[entidx6,1], pch = 20, cex=1.2, col = 'black')
  
  legend("topleft", bg="white", legend=c(entidx1, entidx2, entidx3, entidx4, entidx5, entidx6), 
         lwd = 2, col=c('blue', 'darkgreen', 'red', 'cyan', 'magenta', 'black'))
  
  dev.off()
}  