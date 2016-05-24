pullup <- function(mati, plotting=FALSE, channel=4){
  
  ### first pull up negative peaks
  #for(f in 1:(dim(mati)[2] -1)){
  #s1 <- (mati[,f]* -1)
  #starting <- median(mati[,f])
  #mati[which(s1 > (starting+200))] <- abs(mati[which(s1 > (starting+200))]) * 2
  #}
  ### now the real pullup
  all.cols <- apply(mati[,-channel],1,sum)
  #plot(all.cols, type="l", xlim=c(0,5000))
  #head(mati)
  rows <- dim(all.cols)[1]
  peak.all.cols <- big.peaks.col(all.cols, tre=300)
  #### which are the regions of the peaks to apply pull up, set difference to 100
  vv1 <- which(diff(peak.all.cols$pos) > 100); vv1 <- c(vv1, length(peak.all.cols$pos)) # where peaks cut
  vv2 <- vv1+1; vv2 <- vv2[-length(vv2)]; vv2 <- c(1,vv2)
  regions <- apply(data.frame(cbind(vv2,vv1)), 1, function(x,pos){y <- pos[x[1]:x[2]]; return(y)}, peak.all.cols$pos)
  #########################
  for(i in 1:length(regions)){
    # find which is the maximum for that region in all the channels
    maxis <- apply(mati[,-channel], 2, function(x, regi){max(x[regi])}, regi=regions[[i]])
    ch <- c(1:length(maxis))
    wiii <- which(maxis == max(maxis)) # channel where this peak is real
    rest <- ch[-wiii] # channels where the peak is not real
    # adjust by decreasing 25% the intensity of the real channel
    decrease <- regions[[i]][1]:(regions[[i]][length(regions[[i]])])
    mati[decrease,rest] <- mati[decrease,rest] - (mati[decrease,wiii]*.3)
  }
  
  #all.inds.mats2 <- apply(mati, 2, function(x){x[which(x < 0)] <- 0; return(x)}) 
  all.inds.mats2 <- apply(mati, 2, function(x){x[which(x < 0)] <- x[which(x < 0)]; return(x)}) # THIS PART SEEMS TO BE IMPORTANT FOR DONT GETTING DOUBLE LINES
  #all.inds.mats2  <- data.frame(all.inds.mats2)
  if(plotting == TRUE){
    layout(matrix(1:2,2,1))
    plot(mati[,1], type="l", col="blue")
    colorsh <- c("cornflowerblue", "chartreuse4", "gold2", "red", "orange", "purple")
    for(i in 2:dim(mati)[2]){
      lines(mati[,i], col=colorsh[i])  
    }
    plot(all.inds.mats2[,1], type="l", col="blue")
    colorsh <- c("cornflowerblue", "chartreuse4", "gold2", "red", "orange", "purple")
    for(i in 2:dim(all.inds.mats2)[2]){
      lines(all.inds.mats2[,i], col=colorsh[i])  
    }
  }
  layout(matrix(1,1,1))

  return(all.inds.mats2)
}