`personItemPlot` <- 
  function(raschResult, nBreaks=15, plotTitle="Person Item Histogram",
           xlab = "Relative Frequency", ylab = "Ability",
           col = c("darkgrey","lightgrey"), colTheme, makeLegend=TRUE, 
           legendLabels=c("items", "people"), legendLoc="bottomleft",
           gDevice, file){
  
  if(! missing(colTheme)){
    if(colTheme=="dukes"){
      col <- c("#450084","#CBB677")
    } else if(colTheme=="spartans"){
      col <- c("#003366","#FFCC00")
    } else if(colTheme=="cavaliers"){
      col <- c("#FF7003","#0D3268")
    } else{
      col <- c("darkgrey","lightgrey")
    }	  
  }  
  
  measureRange <- function(items, people){
    mRange <- range(c(items, people))
    mRange[1] <- floor(mRange[1]*10)/10
    mRange[2] <- ceiling(mRange[2]*10)/10
    mRange
  }
  
  mRange <- measureRange(raschResult$item.par$delta.i,raschResult$person.par$theta)
  breakArray <- seq(mRange[1], mRange[2], (mRange[2]-mRange[1])/nBreaks)
  
  itemCounts <- hist(raschResult$item.par$delta.i,plot=FALSE,breaks=breakArray)
  personCounts <- hist(raschResult$person.par$theta,plot=FALSE,breaks=breakArray)
  
  itemCounts$density <- itemCounts$density/sum(itemCounts$density)
  personCounts$density <- personCounts$density/sum(personCounts$density)
  
  maxPeople <- max(personCounts$density)
  maxItems <-  max(itemCounts$density)
  
  startScale <- personCounts$mids[1]
  unitInc <- personCounts$mids[2] - personCounts$mids[1]
  
  if(missing(gDevice)) gDevice <- "screen"
  if(gDevice == "screen"){
    dev.new()    
  } else if(gDevice=="jpg" | gDevice=="jpeg"){
    if(missing(file)) file <- "personItemPlot.jpg"
    jpeg(filename=file)
  } else if(gDevice=="png"){
    if(missing(file)) file <- "personItemPlot.png"
    png(filename=file)
  }
  
  plot(1, type="n",xlim=c(-maxItems,maxPeople),ylim=c(0,nBreaks), main=plotTitle,
       xlab=xlab, ylab=ylab, yaxt="n", xaxt="n")
  
  axYPoints <- floor(mRange[1]):floor(mRange[2]) 
  axis(2,axYPoints,at=((1 + (axYPoints - startScale)/unitInc)))
  
  axXPoints <- seq(ceiling(maxItems*10)/10,0,by=-.1)
  atXPoints <- seq(-ceiling(maxItems*10)/10,0,by=.1) 
  axis(1,axXPoints,at=atXPoints)
  
  axXPoints <- seq(.1, ceiling(maxPeople*10)/10,by=.1)
  atXPoints <- seq(.1, ceiling(maxPeople*10)/10,by=.1) 
  axis(1,axXPoints,at=atXPoints)     
  
  barplot(-itemCounts$density, col=col[1] , horiz=TRUE, space=0, add=TRUE, axes=FALSE)
  barplot(personCounts$density, col=col[2], horiz=TRUE, space=0, add=TRUE, axes=FALSE)
  
  if(makeLegend){
    legend(legendLoc, legendLabels, col=col, pch=15, cex=1.2)
  }
  if(gDevice != "screen") dev.off()
  
}
