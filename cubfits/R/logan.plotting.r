plotCUB.NSE <- function(reu13.df.obs, bMat=NULL, bVec=NULL, phi.bin, n.use.samples=2000,
                     main="CUB", model.label=c("True Model"), model.lty=1,
                     delta_a12=0, a_2=1, positions=c(100, 200, 400, 800), weightedCenters=TRUE, logBins)
{
  ### Arrange data.
  aa.names <- names(reu13.df.obs)
  phi.bin.lim <- range(phi.bin)#range(c(phi.bin, phiMat))
  if(is.null(bMat))
  {
    Eb <- bVec  
  }else{
    lbound <- max(0, length(bMat)-n.use.samples)
    ubound <- length(bMat)
    b.mat <- do.call(cbind, bMat[lbound:ubound]) 
    Eb <- rowMeans(b.mat)  
  }
  Eb <- convert.bVec.to.b(Eb, aa.names)
  
  ### Compute.
  ret.phi.bin <- prop.bin.roc(reu13.df.obs, phi.bin,weightedCenters=weightedCenters, logBins=logBins)
  prediction <- prop.model.nse(Eb, reu13.df.obs, phi.bin.lim, delta_a12=delta_a12, a_2=a_2, positions=positions)
  
  
  ### Fix xlim at log10 scale. 
  lim.bin <- range(log10(ret.phi.bin[[1]]$center))
  lim.model <- range(log10(prediction[[1]][[1]]$center))
  xlim <- c(lim.bin[1] - (lim.bin[2] - lim.bin[1]) / 4,
            max(lim.bin[2], lim.model[2]))
  
  
  mat <- matrix(c(rep(1, 4), 2:21, rep(22, 4)),
                nrow = 7, ncol = 4, byrow = TRUE)
  mat <- cbind(rep(23, 7), mat, rep(24, 7))
  nf <- layout(mat, c(3, rep(8, 4), 2), c(3, 8, 8, 8, 8, 8, 3), respect = FALSE)
  ### Plot title.
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.6, main)
  text(0.5, 0.4, date(), cex = 0.6)
  
  ### Plot results.
  for(i.aa in 1:length(aa.names))
  {
    tmp.obs <- ret.phi.bin[[i.aa]]
    tmp.nse <- lapply(1:length(prediction), function(i.pos){  prediction[[i.pos]][[i.aa]] })
    
    plotbin.NSE(tmp.obs, tmp.nse, main = "", xlab = "", ylab = "",
            lty = model.lty, axes = FALSE, xlim = xlim)
    box()
    main.aa <- oneLetterAAtoThreeLetterAA(aa.names[i.aa])
    text(0, 1, main.aa, cex = 1.5)
    if(i.aa %in% c(1, 5, 9, 13, 17)){
      axis(2)
    }
    if(i.aa %in% 16:19){
      axis(1)
    }
    if(i.aa %in% 1:4){
      axis(3)
    }
    if(i.aa %in% c(4, 8, 12,16)){
      axis(4)
    }
    axis(1, tck = 0.02, labels = FALSE)
    axis(2, tck = 0.02, labels = FALSE)
    axis(3, tck = 0.02, labels = FALSE)
    axis(4, tck = 0.02, labels = FALSE)
  }
  
  ## adding a histogram of phi values to plot
  hist.values <- hist(log10(phi.bin), plot=FALSE, nclass=30)
  plot(hist.values, axes = FALSE, main="", xlab = "", ylab = "")
  legend(x="right", legend=positions, lty = c(1, 2, 4, 3), title="Position")
  box()
  axis(1)

  ### Add label.
#  plot(NULL, NULL, axes = FALSE, main = "", xlab = "", ylab = "",
#        xlim = c(0, 1), ylim = c(0, 1))
#  legend(0.1, 0.8, model.label, lty = model.lty, box.lty = 0)
  
  ### Plot xlab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Production Rate (log10)")
  
  ### Plot ylab.
  plot(NULL, NULL, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE)
  text(0.5, 0.5, "Propotion", srt = 90)
}


### This plots one amino acid the NSE model.
plotbin.NSE <- function(ret.bin, ret.model = NULL, main = NULL,
    xlab = "Production Rate (log10)", ylab = "Proportion",
    xlim = NULL, lty = 1, x.log10 = TRUE, stderr = FALSE, ...){
  #browser();
  if(x.log10){
    ret.bin$center <- log10(ret.bin$center)
    if(!is.null(ret.model)){
      for(ii in 1:length(ret.model)){ ret.model[[ii]]$center <- log10(ret.model[[ii]]$center) }
    }
  }

  ### Observed dots and whiskers.
  if(is.null(xlim)){
    if(!is.null(ret.model)){
      # x.lim <- range(c(ret.bin$center, ret.model$center))
      lim.bin <- range(ret.bin$center)
      lim.model <- range(ret.model[[1]]$center)
      x.lim <- c((lim.bin[1] + lim.model[1]) / 2,
                 max(lim.bin[2], lim.model[2]))
    } else{
      x.lim <- range(ret.bin$center)
    }
  } else{
    x.lim <- xlim
  }
  y.lim <- c(0, 1) + c(-0.05, 0.05)

  u.codon <- sort(unique(as.character(ret.bin$codon)))
  u.center <- unique(ret.bin$center)
  color <- get.color(u.codon)
  color.alpha <- get.color(u.codon, color = .CF.PT$color.alpha)
  ### Reorder R for better legend.
  if(all(u.codon %in% .CF.GV$synonymous.codon$R)){
    u.codon <- u.codon[c(3:6, 1:2)]
    color <- color[c(3:6, 1:2)]
    color.alpha <- color.alpha[c(3:6, 1:2)]
  }

  ngenes.codon.total <- 0
  for(i.codon in 1:length(u.codon)){  
    ngenes.codon <- ret.bin[ret.bin$codon == u.codon[i.codon], ]
    ngenes.codon.total <- ngenes.codon.total + ngenes.codon$codon.count
  }
  ret.bin.codon <- ret.bin[ret.bin$codon == u.codon[1],]
  plot(ret.bin.codon$center, ngenes.codon.total, ,type="h", pch=22, col="#D8D8D888", lwd=7, xaxt="n", yaxt="n", bty="n", xlab="",ylab="", xlim=x.lim)    
  par(new=TRUE)
  
  
  ### Make an empty plot
  plot(NULL, NULL, xlim = x.lim, ylim = y.lim,
       main = main, xlab = xlab, ylab = ylab, ...)

  ### Add observed dots for means and whiskers for standard deviations.
  for(i.codon in 1:length(u.codon)){
    ret.bin.codon <- ret.bin[ret.bin$codon == u.codon[i.codon],]

    ### Add vertical lines.
    for(i.center in 1:nrow(ret.bin.codon)){
      if(!is.na(ret.bin.codon$freq.std[i.center])){
        if(stderr){
          freq.bar <- ret.bin.codon$freq.stderr[i.center]
        } else{
          freq.bar <- ret.bin.codon$freq.std[i.center]
        }
        tmp.y <- ret.bin.codon$freq.mean[i.center] + c(-1, 1) * freq.bar
        tmp.y[tmp.y < 0] <- 0
        tmp.y[tmp.y > 1] <- 1
        lines(
          list(x = rep(ret.bin.codon$center[i.center], 2),
               y = tmp.y),
          col = color.alpha[i.codon], lwd = 0.8)
        
      }
    }
    ### Add points.
    points(ret.bin.codon$center,
           ret.bin.codon$freq.mean,
           pch = 19, col = color[i.codon], cex = 0.5)    
  }
  
  
  ### Add modeled lines.
  u.codon.star <- u.codon
  if(!is.null(ret.model)){
    plotaddmodel(ret.model, lty, u.codon, color, x.log10 = FALSE)

    ### Add focal codon.
    u.codon.star <- attr(ret.model[[1]], "u.codon.star")
    if(!is.null(u.codon.star)){
      if(all(u.codon %in% .CF.GV$synonymous.codon$R)){
        u.codon.star <- u.codon.star[c(3:6, 1:2)]
      }
    }
  }

  ### Add legends.
  legend(x.lim[1], y.lim[2], u.codon.star, col = color, bg = NULL,
         box.lty = 0, lty = 1, pch = 19, cex = 0.8)
} # End of plotbin().

