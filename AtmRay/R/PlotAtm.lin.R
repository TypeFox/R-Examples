PlotAtm.lin = function(ATM, zlim = c(0, 100), winddir = 90, col = sky.colors(500), TOPO = NULL){
  ###############################
  # define scale image function
  
  scale_image=function (x = seq(from = 0, to = 1, length.out = dim(Z)[1]), 
    y = seq(from = 0, to = 1, length.out = dim(Z)[2]), Z, l = 0.5, 
    w = 0.075, title = "", toss.outliers = 2, log = FALSE, 
    center = NULL, pretty = 0, rescale=TRUE, correction=1, mar=NULL, breaks=NULL,col=heat.colors(60),scale.lab=NULL,...){
    if(is.null(mar)){par(mar = c(5.1, 4.1, 4.1, 4.1))}
    if(rescale){#zlim=correction*zlim/mean(Z[is.finite(Z)])
      Z=correction*Z/mean(Z[is.finite(Z)],na.rm=TRUE);center=1}
    
    if (log) {
      Z = log(Z, 10)
      if(!is.null(breaks)){breaks = log(breaks, 10)}
    }
    
    if(!is.matrix(Z) & (length(Z) == (length(x)-1) * (length(y)-1)) & length(x) & length(y)){
      Z = matrix(Z,length(x)-1,length(y)-1)
    }
    
    Z[Z==Inf]=max(Z[is.finite(Z)],na.rm=TRUE)
    Z[Z==-Inf]=min(Z[is.finite(Z)],na.rm=TRUE)
    
    
    if (toss.outliers) {
      Z[Z < (mean(Z, na.rm = TRUE) - toss.outliers * sd(as.vector(Z), na.rm = TRUE))] = mean(Z, 
         na.rm = TRUE) - toss.outliers * sd(as.vector(Z), na.rm = TRUE)
      Z[Z > (mean(Z, na.rm = TRUE) + toss.outliers * sd(as.vector(Z), na.rm = TRUE))] = mean(Z, 
         na.rm = TRUE) + toss.outliers * sd(as.vector(Z), na.rm = TRUE)
    }
    
    m = dim(Z)[1]
    n = dim(Z)[2]
    L = ceiling(n * l)
    H = w/(1 - w)
    
    if(is.null(breaks)){
      breaks = seq(from = min(Z[is.finite(Z)], na.rm = TRUE), to = max(Z[is.finite(Z)], na.rm = TRUE), length.out = 1 + length(col))
    }else if(length(breaks)!=length(col)+1){stop('length(breaks) != length(col)+1')}
    scale = approx(1:length(breaks)/length(breaks),breaks,1:L/L)$y
    S = rbind(rep(NA, n), c(rep(NA, floor((n - L)/2)), scale, 
      rep(NA, ceiling((n - L)/2))), rep(NA, n))
    Z = rbind(Z, S)
    x = c(x, max(x) + 5 * (x[2] - x[1]), max(x) + H * (max(x) - 
      x[1]), max(x) + H * (max(x) - x[1]) + 3 * (x[2] - x[1]))
    
    image(x, y, Z, breaks=c(-10^10,breaks[2:(length(breaks)-1)],10^10),col=col,...)
    
    if(is.null(scale.lab)){
      if (log) {
        scale.lab = signif(10^c(min(breaks), max(breaks)), 3)
      }
      if (!log) {
        scale.lab = signif(c(min(breaks), max(breaks)), 3)
      }
      scale.lab = c(scale.lab[1], center, scale.lab[2])
      
      
      if (log & !is.null(center)) {
        center = log(center, 10)
      }
      if (pretty & !log) {
        scale.lab = pretty(scale.lab, pretty)
      }
      if (pretty & log) {
        scale.lab = signif(10^pretty(log(scale.lab, 10), pretty), 2)
      }
      print(scale.lab)
    }
    if (!log) {
      at = y[as.vector(sapply(scale.lab, function(x) {
        which.min(abs(Z[m + 2, ] - x))
      }))]
    }
    if (log) {
      at = y[as.vector(sapply(scale.lab, function(x) {
        which.min(abs(Z[m + 2, ] - log(x, 10)))
      }))]
    }
    print(at)
    axis(side = 4, at = at, labels = scale.lab, las = 0)
    mtext(title, side = 4, line = 2.25, cex = par("cex"))
  }
  ################################################
  # done defining scale image
  
  scale_image(x = -100:100, y = seq(min(zlim), max(zlim), length.out = 501), Z = matrix(seq(min(zlim), max(zlim), length.out = 500) - ATM$z0, 200, 500, byrow = TRUE) * ATM$gc + ATM$c0, col = col, xaxt = 'n', xlab = '', ylab = 'Elevation (m)', rescale = FALSE, l = 0.75, title = 'Intrinsic Sound Speed (m/s)')
  
  w0 = ATM$wx0 * sin(winddir * pi/180) + ATM$wy0 * cos(winddir * pi/180)
  gw = ATM$gwx * sin(winddir * pi/180) + ATM$gwy * cos(winddir * pi/180)
  
  z = seq(min(zlim), max(zlim), length.out = 7)[2:6]
  w = w0 + gw * (z - ATM$z0)
  
  if(sum(w != 0) > 0){
    arrows(0, z, w/40 * 100, z, length = 0.125, lwd = 2)
    text(w/40 * 100, z, paste(round(w, 1), 'm/s'), adj = (1 - sign(w))/2)
  }
  
  if(!is.null(TOPO)){
    polygon(c(-100, RESCALE(TOPO$e, -100, 100, min(TOPO$e), max(TOPO$e)), 100), c(0, TOPO$z[,which.min(abs(TOPO$n))], 0), 100, 90, col = 'darkblue', border = 'black')
  }
}
