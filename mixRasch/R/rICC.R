`rICC` <- function(delt, theta, itemVector, xlim, ylim, plotTitle, 
                   xlab, ylab, col = c("black","white"), colTheme, 
                   expectedScore=FALSE, empICC=FALSE, empOnly=FALSE, gDevice, 
                   file, ...){
  
  
  cutFind <- function(xxx, nCuts){
    
    rng <- range(xxx, na.rm=TRUE)
    
    breakCheck <- function(){
      cutPlus <- (rng[2] - rng[1])/nCuts
      breaks <- rng[1] + c(0,cumsum(rep(cutPlus,nCuts)))
      mids <- 1:nCuts
      for(j in 1:nCuts){
        mids[j] <- (breaks[j] + breaks[j+1])/2
      }
      breaks[1] <- breaks[1] - .1
      breaks[nCuts+1] <- breaks[nCuts+1] + .2
      cats <- cut(xxx,breaks)
      list(cats, mids, sort(table(cats))[1:2])
    }
    if(! missing(nCuts)){
      out <- breakCheck()
    }else{
      for(nCuts in 3:20){
        res <- breakCheck()
        out <- res[[1]]
        if(res[[3]][1] < 15 | res[[3]][2] < 15){ 
          nCuts <- nCuts - 1
          out <- breakCheck()
          break
        }
      }  	
    }
    out  
  }		
  if(! missing(colTheme)){
    if(colTheme=="dukes"){
      col <- c("#450084","#CBB677")
    } else if(colTheme=="spartans"){
      col <- c("#003366","#FFCC00")
    } else if(colTheme=="cavaliers"){
      col <- c("#0D3268","#FF7003")
    } else if(colTheme=="greys"){
      col <- c("black","lightgrey")
    } else {
      col <- c("black","white")	
    }	  
  }
  
  if(missing(ylim)) ylim <- c(0,1)
  
  if(! missing(empOnly)){
    if(empOnly) empICC <- TRUE
  }
  
  if(missing(plotTitle)) plotTitle <- "Item Characteristic Curve"
  if(missing(xlab)) xlab <- "Ability"
  if(missing(ylab)) ylab <- "Probability"
  
  nCat <- length(delt) + 1
  
  if(nCat > 2 & empICC) expectedScore <- TRUE
  
  if(missing(theta)){
    if(empICC){
      warning("You did not provide estimated Theta values. Your empirical ICC will not be created.")
      empICC <- FALSE
      empOnly <- FALSE
    }
    if(missing(xlim)){
      xlim <- c(-3,3)
    }		
  } else{
    if(missing(xlim)) xlim <- range(theta)
  }
  
  if(empICC){
    if(missing(itemVector)){
      warning("You did not provide the observed item vector. Your empirical ICC will not be created.")
      empICC <- FALSE
      empOnly <- FALSE
    } else{
      theCut <- cutFind(theta)
      itemVector <- na.omit(cbind(itemVector,theCut[[1]]))
      itemMeans <- c(by(itemVector[,1],itemVector[,2],mean,na.rm=TRUE))
      meansLoc <- theCut[[2]]
    }
  }
  
  if(missing(gDevice)) gDevice <- "screen"
  if(gDevice == "screen"){
    #dev.new()    
  } else if(gDevice=="jpg" | gDevice=="jpeg"){
    if(missing(file)) file <- "iccPlot.jpg"
    jpeg(filename=file)
  } else if(gDevice=="png"){
    if(missing(file)) file <- "iccPlot.png"
    png(filename=file)
  }
  
  if(! empOnly){	
    values <- seq(xlim[1],xlim[2],.05)
    probs <- P.xj(delt, values)
    
    if(nCat > 2){
      probs <- rbind(1-colSums(probs), probs)
      if(expectedScore){
        #if(missing(ylim)) 
        ylim <- c(0, nCat-1)
        probs <- matrix(colSums(probs*seq(0,nCat-1,1)), nrow=1)
        if(missing(ylab)) ylab <- "Expected Item Score"    	
      }   	
    }
    
    plot(values, probs[1,], col=col[1], type="l", ylim=ylim, main=plotTitle, xlab=xlab, ylab=ylab, ...)
    
    if(nCat > 2 & ! expectedScore){
      for(i in 2:nCat){
        points(values, probs[i,], col=col[1], type="l", ...)
      }
    }
    if(empICC){
      points(meansLoc, itemMeans, type="p", col=col[2], pch=16, ...)
      points(meansLoc, itemMeans, type="b", col=col[1], lty=2, ...)
    }
  } else{
    plot(meansLoc, itemMeans, type="p", col=col[2], pch=16, ylim=c(0,1), main=plotTitle, xlab=xlab, ylab=ylab, ...)
    points(meansLoc, itemMeans, type="b", col=col[1], lty=2, ...)
  }	 	
  if(gDevice != "screen") dev.off()
}
