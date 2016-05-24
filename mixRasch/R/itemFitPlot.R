`itemFitPlot` <- function(raschResult, itemSet, useItemNames = TRUE, 
                        fitStat = "infit", plotTitle="Item Fit Plot",
                        xlab, ylab, xlim, ylim, refLines, 
                        col = c("black","white"), colTheme,
                        gDevice, file){
  
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
  
  if(missing(ylab)) ylab <- "Difficulty" 
  
  if(missing(itemSet)){
    itemSet <- 1:raschResult$i.stat$n.i
  } else{
    nameSet <- colnames(raschResult$item.par$delta) %in% itemSet
    if(any(nameSet)) itemSet <- nameSet
  }	
  
  fitV <- raschResult$item.par$in.out[itemSet,fitStat]
  
  if(fitStat=="infit"){
    refBump <- .1
    if(missing(xlab)) xlab <- "Infit"
    if(missing(refLines)) refLines <- c(.7,1.3)
    #if(missing(xlim)) xlim <- c(.3,1.7)
  } else if(fitStat=="outfit"){
    refBump <- .1
    if(missing(xlab)) xlab <- "Outfit"
    if(missing(refLines)) refLines <- c(.7,1.3)
    #if(missing(xlim)) xlim <- c(.3,1.7)
  } else if(fitStat=="in.Z"){
    refBump <- .2
    if(missing(xlab)) xlab <- "Standardized Infit"
    if(missing(refLines)) refLines <- c(-2,2)
    #if(missing(xlim)) xlim <- c(-3,3)  	
  } else{
    refBump <- .2
    fitStat <- "out.Z"
    if(missing(xlab)) xlab <- "Standardized Outfit"
    if(missing(refLines)) refLines <- c(-2,2)
    #if(missing(xlim)) xlim <- c(-3,3)
  }	
  
  if(missing(xlim)){
    xlim <- range(fitV)
    xlim[1] <- xlim[1] - .1
    xlim[2] <- xlim[2] + .1
    if(xlim[1] > refLines[1]) xlim[1] <- refLines[1] - refBump
    if(xlim[2] < refLines[2]) xlim[2] <- refLines[2] + refBump	
  } 
  
  if(useItemNames){
    itemNames <- colnames(raschResult$item.par$delta)[itemSet]
  } else{
    itemNames <- paste(itemSet)
  }	  
  
  myWidth <- 4
  if(missing(gDevice)) gDevice <- "screen"
  if(gDevice == "screen"){
    dev.new(width=myWidth)    
  } else if(gDevice=="jpg" | gDevice=="jpeg"){
    if(missing(file)) file <- "itemFitPlot.jpg"
    jpeg(width=480*myWidth/5, filename=file)
  } else if(gDevice=="png"){
    if(missing(file)) file <- "itemFitPlot.png"
    png(width=480*myWidth/5, filename=file)
  }
  plot(fitV, raschResult$item.par$delta.i[itemSet],
       xlim=xlim, type="n", xlab=xlab, ylab=ylab, main=plotTitle)  
  polygon(c(refLines[2],refLines[2],500,500),c(-500,500,500,-500),col=col[2])
  polygon(c(refLines[1],refLines[1],-500,-500),c(-500,500,500,-500),col=col[2])  
  text(fitV, raschResult$item.par$delta.i[itemSet], col=col[1], labels=itemNames)
  abline(v=refLines)
  if(gDevice != "screen") dev.off()
}