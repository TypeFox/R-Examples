`cttICC` <- function(scores, itemVector, xlim, ylim, plotTitle, 
                   xlab, ylab, col = c("black","white"), colTheme,
				   gDevice, file, ...){
  
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
    
  if(missing(plotTitle)) plotTitle <- "Item Characteristic Curve"
  if(missing(xlab)) xlab <- "Test Score"
  if(missing(ylab)) ylab <- "Item Mean"
    
  if(missing(xlim)) xlim <- range(scores,na.rm=TRUE)
  if(missing(ylim)) ylim <- c(min(itemVector,na.rm=TRUE), max(itemVector,na.rm=TRUE)) 
  theCut <- cutFind(scores)
  itemVector <- na.omit(cbind(itemVector,theCut[[1]]))
  itemMeans <- c(by(itemVector[,1],itemVector[,2],mean,na.rm=TRUE))
  meansLoc <- theCut[[2]]
  
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
  
  plot(meansLoc, itemMeans, type="p", col=col[2], pch=16, ylim=ylim, main=plotTitle, xlab=xlab, ylab=ylab, ...)
  points(meansLoc, itemMeans, type="b", col=col[1], lty=2, ...)
  	 	
  if(gDevice != "screen") dev.off()
}
