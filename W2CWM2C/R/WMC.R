WMC <-
function(inputDATA, Wname, J, device="screen", filename,
          Hfig, WFig, Hpdf, Wpdf) { 
  
 if(is.ts(inputDATA) != "TRUE")
  stop("The input data is not a time series, please check the ts 
  function in the R manual pages. Bye, thank you for your interest 
  in our program. \n")

 NAMES <- colnames(inputDATA) 

 #: Compute the dimensions 
  MNL <- dim(inputDATA)
  M  <- MNL[2] #No. columns
  N  <- MNL[1] #No. elements
  if(M >= N) stop("Be careful with the input data, there 
   are more columns (variables) than number of elements.")

  #Comp. the MODWT
  LIST <- as.list("NULL")    
  for (l in 1:M) { 
    tmp.wt  <- modwt(inputDATA[,l], Wname, J) 
    tmp.bw  <- brick.wall(tmp.wt, Wname)	
    LIST[l] <- list(tmp.bw)
  }

  LS                <- wave.multiple.correlation(LIST, N) 
  returns.modwt.cor <- LS$xy.mulcor[1:J,]
  YmaxR             <- LS$YmaxR

   ## Devices options: png & jpg; esp & pdf! 
  if (device=="png") {
   fileout <- paste("WMC_", filename, ".png", sep="")
   png(fileout, height=Hfig, width=WFig)
  }

  if (device=="jpeg" || device=="jpg") {
   fileout <- paste("WMC_", filename, ".jpg", sep="")
   jpeg(fileout, height=Hfig, width=WFig)
  }

  if (device=="pdf") {
   fileout <- paste("WMC_", filename, ".pdf", sep="")
   pdf(fileout, height=Hpdf, width=Wpdf)
  }

  if (device=="eps") {
   fileout <- paste("WMC_", filename, ".eps", sep="")
   postscript(fileout, height=Hpdf, width=Wpdf)
  }

  par(mfrow=c(1,1), las=0, mar=c(5,4,4,2)+.1)
  matplot(2^(0:(J-1)), returns.modwt.cor[-(J+1),], type="b",
  log="x", pch="*LU", xaxt="n", lty=1, col=c(1,4,4),
  xlab="Wavelet Scale", ylab="Wavelet Multiple Correlation", 
  ylim=c(-0.1,1))
  axis(side=1, at=2^(0:7))
  abline(h=0)
  text(2^(0:7),rep(-0.1,J),labels=NAMES[YmaxR[-1]],adj=0.5,cex=.85)

  if (device != "screen")
  dev.off()
  
  cat("labels",  NAMES[YmaxR[-1]])
 
  return(list(LS=LS))

}
