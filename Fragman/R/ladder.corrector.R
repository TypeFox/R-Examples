ladder.corrector <- function(stored, to.correct, ladder, thresh=200, env = parent.frame()){
  list.data <- env$list.data.covarrubias
  
  vvv <- which(names(list.data) %in% to.correct)
  chan <- dim(stored[[1]])[2]
  ###
  #--
  ###
  for(j in 1:length(vvv)){
    cat(paste("\nSTARTING SAMPLE",j,"\n"))
    www <- vvv[j]
    x <- stored[[www]][,chan]
    roxy <- big.peaks.col(x, thresh)
    ylimi <- 2 * median(sort(roxy$hei, decreasing=TRUE)) 
    plot(x, type="l", ylim=c(0,ylimi), main=to.correct[j])
    cat(paste("\nYou need you to click over peaks corresponding to the following ladder: \n\n", paste(ladder, collapse = " "), "\n\nPress 'esc' when you are done. \n"))
    good <- locator(type="p", pch=20, col="red")$x
    #############################################
    ## now use such values to correct the model
    corrected <- numeric()
    for(k in 1:length(good)){
      yoyo <- abs(roxy$pos - good[k])
      corrected[k] <- roxy$pos[which(yoyo == min(yoyo))]
    }
    corrected.hei <- roxy$hei[which(roxy$pos %in% corrected)]
    if(length(corrected) != length(ladder)){
      cat("ERROR!! you did not select the same number of peaks than the number of elements of your ladder ")
    }else{
      list.data[[www]]$pos <- corrected
      list.data[[www]]$hei <- corrected.hei
      list.data[[www]]$wei <- ladder
      list.data[[www]]$corr <- cor(ladder, corrected)
    }
    ##############################################
  }
  ###
  #--
  ###
  for(l in 1:length(vvv)){
    www <- vvv[l]
    x <- stored[[www]][,chan]
    limi <- sort(list.data[[www]]$hei, decreasing = TRUE)
    plot(x, type="l", xaxt="n", ylim=c(0,(limi[3]+1000)), cex.axis=0.6, las=2, xlim=c((min(list.data[[www]]$pos)-100),(max(list.data[[www]]$pos)+100)), col=transp("grey35",0.7), ylab="RFU", xlab="", lwd=2)
    axis(1, at=list.data[[www]]$pos, labels=list.data[[www]]$wei, cex.axis=0.6)
    points(x=list.data[[www]]$pos, y=list.data[[www]]$hei,cex=1.1, col=transp("black",0.85))
    points(x=list.data[[www]]$pos, y=list.data[[www]]$hei, pch=20, col=transp("red",0.7))
    legend("topleft", legend=paste("Correlation:",round(list.data[[www]]$corr, digits=4), sep=""), bty="n")
    legend("topright", legend=c("Peaks selected"), col=c("red"), bty = "n", pch=c(20), cex=0.85)
    
  }
  env$list.data.covarrubias <- list.data
  cat("\nJOB DONE!!! ladder has been corrected for your samples")
}