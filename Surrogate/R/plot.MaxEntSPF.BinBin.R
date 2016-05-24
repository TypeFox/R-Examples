plot.MaxEntSPF.BinBin <- function(x, SPF.Fit, Type="All.Histograms", Col="grey", ...){
                            
  Object <- x 

  if (missing(Col)) {Col = "grey"}
  
  if (Type=="All.Histograms"){
    
    plot(0:100, 0:100, axes=F, xlab="", ylab="", type="n", xlim=c(0, 1), ...)  
    par(mfrow=c(3, 3), mar = c(4.5, 7, 4, 1), oma=rep(0, times=4))  
    
    # T = -1, S = -1
    if (length(unique(SPF.Fit$r_min1_min1)) > 1){
      hist(SPF.Fit$r_min1_min1, main="", col=Col, xlim=c(0, 1),
           xlab=expression(r(-1,-1)), cex.lab=1.3)}
    if (length(unique(SPF.Fit$r_min1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_min1_min1, lwd=2)
    
    # T = 0, S = -1
    if (length(unique(SPF.Fit$r_0_min1)) > 1){
      hist(SPF.Fit$r_0_min1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,-1)))}
    if (length(unique(SPF.Fit$r_0_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_0_min1, lwd=2)
    
    # T = 1, S = -1
    if (length(unique(SPF.Fit$r_1_min1)) > 1){
      hist(SPF.Fit$r_1_min1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,-1)))}
    if (length(unique(SPF.Fit$r_1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_1_min1, lwd=2)
    
    # T = -1, S = 0
    if (length(unique(SPF.Fit$r_min1_0)) > 1){
      hist(SPF.Fit$r_min1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(-1,0)))}
    if (length(unique(SPF.Fit$r_min1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_min1_0, lwd=2)
    
    # T = 0, S = 0
    if (length(unique(SPF.Fit$r_0_0)) > 1){
      hist(SPF.Fit$r_0_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,0)))}
    if (length(unique(SPF.Fit$r_0_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_0_0, lwd=2)
    
    # T = 1, S = 0
    if (length(unique(SPF.Fit$r_1_0)) > 1){
      hist(SPF.Fit$r_1_0, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,0)))}
    if (length(unique(SPF.Fit$r_1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_1_0, lwd=2)
    
    # T = -1, S = 1
    if (length(unique(SPF.Fit$r_min1_1)) > 1){
      hist(SPF.Fit$r_min1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(-1,1)))}
    if (length(unique(SPF.Fit$r_min1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_min1_1, lwd=2)
    
    # T = 0, S = 1
    if (length(unique(SPF.Fit$r_0_1)) > 1){
      hist(SPF.Fit$r_0_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,1)))}
    if (length(unique(SPF.Fit$r_0_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_0_1, lwd=2)
    
    
    # T = 1, S = 1
    if (length(unique(SPF.Fit$r_1_1)) > 1){
      hist(SPF.Fit$r_1_1, main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,1)))}
    if (length(unique(SPF.Fit$r_1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_1_1, lwd=2)
            }
    
    

  if (Type=="All.Densities"){
    
    plot(0:100, 0:100, axes=F, xlab="", ylab="", type="n", xlim=c(0, 1), ...)  #LS!
    par(mfrow=c(3, 3), mar = c(4.5, 7, 4, 1), oma=rep(0, times=4), xpd=FALSE)  #LS!
    
    # T = -1, S = -1
    if (length(unique(SPF.Fit$r_min1_min1)) > 1){
      plot(density(SPF.Fit$r_min1_min1, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(-1,-1)))
    }
    if (length(unique(SPF.Fit$r_min1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_min1_min1, lwd=2)
    
    
    mtext(side = 3, expression(paste(Delta, "T = -1")), cex=1.5, padj = -1.6)   
    mtext(side = 2, expression(paste(Delta, "S = -1")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = -1
    if (length(unique(SPF.Fit$r_0_min1)) > 1){
      plot(density(SPF.Fit$r_0_min1, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,-1)))}
    if (length(unique(SPF.Fit$r_0_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    abline(v = Object$r_0_min1, lwd=2)
    mtext(side = 3, expression(paste(Delta, "T = 0")), cex=1.5, padj = -1.6)   
    
    # T = 1, S = -1
    if (length(unique(SPF.Fit$r_1_min1)) > 1){
      plot(density(SPF.Fit$r_1_min1, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,-1)))}
    if (length(unique(SPF.Fit$r_1_min1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    abline(v = Object$r_1_min1, lwd=2)
    mtext(side = 3, expression(paste(Delta, "T = 1")), cex=1.5, padj = -1.6)   
    
    
    # T = -1, S = 0
    if (length(unique(SPF.Fit$r_min1_0)) > 1){
      plot(density(SPF.Fit$r_min1_0, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(-1,0)))}
    if (length(unique(SPF.Fit$r_min1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    abline(v = Object$r_min1_0, lwd=2)
    mtext(side = 2, expression(paste(Delta, "S = 0")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = 0
    if (length(unique(SPF.Fit$r_0_0)) > 1){
      plot(density(SPF.Fit$r_0_0, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,0)))}
    if (length(unique(SPF.Fit$r_0_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_0_0, lwd=2)
    
    
    # T = 1, S = 0
    if (length(unique(SPF.Fit$r_1_0)) > 1){
      plot(density(SPF.Fit$r_1_0, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,0)))}
    if (length(unique(SPF.Fit$r_1_0)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_1_0, lwd=2)
    
    
    # T = -1, S = 1
    if (length(unique(SPF.Fit$r_min1_1)) > 1){
      plot(density(SPF.Fit$r_min1_1, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(-1,1)))}
    if (length(unique(SPF.Fit$r_min1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    
    abline(v = Object$r_min1_1, lwd=2)
    
    mtext(side = 2, expression(paste(Delta, "S = 1")), cex=1.5, padj = -3.6)   
    
    # T = 0, S = 1
    if (length(unique(SPF.Fit$r_0_1)) > 1){
      plot(density(SPF.Fit$r_0_1, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(0,1)))}
    if (length(unique(SPF.Fit$r_0_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_0_1, lwd=2)
    
    
    # T = 1, S = 1
    if (length(unique(SPF.Fit$r_1_1)) > 1){
      plot(density(SPF.Fit$r_1_1, na.rm=T), main=" ", col=Col, cex.lab=1.3, xlim=c(0, 1),
           xlab=expression(r(1,1)))}
    if (length(unique(SPF.Fit$r_1_1)) <= 1){
      plot(x=0, col=0, axes=F, xlab="", ylab= " ")}
    abline(v = Object$r_1_1, lwd=2)
    
    par(mfrow=c(1, 1), c(5, 4, 4, 2) + 0.1)
    
  }
  
}


    