## functions for NY times plots

##################################################
## setup plot window
plotSetup = function(z,
  ylim = range(coredata(z), na.rm=TRUE),
  xlim = range(index(z)),
  bg = gray(.95),                            # background color
  ...) {
  if(!is.null(bg)) {
    par(bg=bg)
  }

  plot.new()
  plot.window(xlim=xlim, ylim=ylim,bg=bg,...)
}
  
## fancyYaxis -- use dashes, turn labels
fancyYaxis = function(z,
  at = pretty(na.omit(as.numeric(z))),
  labels = at,
  las = 2,
  lty=0,                                # for no axes
  padj = -1,                            # put y labels above
  dashes = TRUE,
  dash.col = gray(.6),
  ...
  ) {

  axis(2, at=at, labels = labels, las=las, lty=lty, padj=padj,
       ...)
  if(dashes == TRUE) {
    abline(h=at, lty=2, col=dash.col)
  }
}


## fancyXaxis made with axis()
fancyXaxis = function(z,
  at = pretty(index(z)),
  dateFormat="%Y-%m-%d",
  labels = format(as.Date(at),format=dateFormat),
  lty=1,                                # set to 0 no axes, ticks
  ...
  ) {

  axis(1,at=at, labels=labels,lty=lty,...)
}
##
##
##################################################

### helper functions
sapply.zoo = function(z, f) {
  tmp = sapply(as.data.frame(coredata(z)),f)
  return(zoo(tmp, index(z)))
}

## assume z use POSIXct for index
insertZeroes = function(z) {
  zVals = coredata(z)
  zDates = index(z)                     #assume POSIXct
  tmp = zVals[1]
  tmpDates = zDates[1]
  n = length(z)

  for(i in 1:(n-1)) {
    tmp[length(tmp)+1] = zVals[i]
    tmpDates[length(tmpDates)+1] = zDates[i]
    
    if(zVals[i]*zVals[i+1] < 0) {
      ## insert 0 at intermediate time
      tmp[length(tmp) +1] = 0
      tmpDates[length(tmpDates) +1] =  mean(zDates[i:(i+1)])
    }
  }
  
  tmp[length(tmp)+1] = zVals[n]
  tmpDates[length(tmpDates)+1] = zDates[n]

  return(zoo(tmp,tmpDates))
}

##
##
##################################################

## Make barplot or histogram like
linesBarplot = function(z,
  proportion = .75,
  col = "black",
  border=col,
  ...) {

  xx = index(z)
  yy = coredata(z)
  n = length(z)

  ## recycle color
  n = length(z)
  col = rep(col, length.out= n)
  border = rep(border, length.out= n)
  
  for(i in 1:(n-1)) {
    theDiff = as.numeric(xx[i+1]) - as.numeric(xx[i])
    gap = (1-proportion)/2*theDiff
    rect(xx[i]+gap,0,xx[i+1]-gap,yy[i], col=col[i], border=border[i], ...)
  }

  ## draw last bar using last theDiff, gap
  rect(xx[n]+gap,0,xx[n]+ theDiff - gap, yy[n], col= col[n],
       border=border[n],...)
}

plotBarplot = function(z,
  ...
  ) {
  plotSetup(z)
  linesBarplot(z,...)
  fancyXaxis(z)
  fancyYaxis(z)
}

## how to plot alternate years
plotAlternateYears = function(z,
  oddcol = "gray",
  evencol = "black",
  ...
  ) {
  is.even = function(x) x %/% 2 == x/2
  col = rep(oddcol,length.out = length(z))
  col[is.even(as.numeric(format(index(z),"%Y")))] = evencol

  
  plotBarplot(z,col=col)
}
  
##
##
##################################################

## histogram like with color
linesHistlike = function(z,
  minVal = min(na.omit(coredata(z))),
  col="black",
  border=col,
  ...
  ) {

  xx = index(z)
  yy = coredata(z)
  for(i in 1:(length(z)-1)) {
    rect(xx[i],minVal,xx[i+1],yy[i],
         col=col, border=border,...
         )
  }
}

plotHistlike = function(z,...) {
  plotSetup(z)
  linesHistlike(z,...)
  fancyXaxis(z)
  fancyYaxis(z)
}

##
##
##################################################

## filled in ts plot
## see example in polygon() 
linesFill = function(z,
  minVal = min(coredata(z)),
  col="gray",
  border="black",
  ...) {

  xx = c(index(z),rev(index(z)))
  yy = c(coredata(z),rep(minVal,length=length(coredata(z))))
  polygon(xx, yy, col=col, border=border)

}

plotFill = function(z,...) {
  plotSetup(z)
  linesFill(z,...)
  fancyXaxis(z)
  fancyYaxis(z)
}






## plot each x differently -- OHLC plot
#### real function is in tseries package

## z has Open High Low Close as matrix columns
makeOHLCforOne = function(z, leftx,rightx,...) {
    midx = mean(c(leftx,rightx))
    segments(midx,z[,'Low'], midx, z[,'High'],...)
    segments(leftx,z[,'Open'], midx, z[,'Open'],...)
    segments(midx,z[,'Close'], rightx, z[,'Close'],...)
  }

plotOHLCgraph = function(z,
  lwd = 3,
  ...) {
  ## zoo object with variables Oen, High, Low, Close
  plotSetup(z,...)
  fancyXaxis(z)
  fancyYaxis(z)
  n = nrow(z)

  zDates = index(z)
  
  for(i in 1:(n-1))  {
    makeOHLCforOne(z[i,], zDates[i], zDates[i+1], lwd=lwd)
  }
}



## inset plot
## use par(fig=c(x1,x2,y1,y2), mai=c(0,0,0,0), new=TRUE) to set up inset
## store are restore oldpar settings
insetPlotExample = function(z.ohlc) {

  z = z.ohlc[,'Close']
  
  plotSetup(z)
  lines(z)
  fancyXaxis(z)
  fancyYaxis(z)

  op = par(no.readonly=T)
  on.exit({par(op)})

  par(fig=c(.6,.8,.6,.8), mai=.05+c(0,0,0,0), new=TRUE)
  z.ohlc = z.ohlc[5:10,]                 # fixed for data set
  plot.window(xlim = range(index(z.ohlc)), ylim = range(coredata(z.ohlc), na.rm=TRUE))
  box()
  zDates = index(z.ohlc)
  n = length(zDates)
  for(i in 1:(n-1))  {
    makeOHLCforOne(z.ohlc[i,], zDates[i], zDates[i+1], lwd=3)
  }

  fancyXaxis(at = zDates[1:5],labels=format(as.Date(zDates[1:5]),"%a"))

}


###### show red and black

plotAboveBelow = function(z,
  acol = "black",
  bcol = "red",
  ...
  ) {

  plotSetup(z)
  fancyXaxis(z)
  fancyYaxis(z)

  znew = insertZeroes(z)

  tmp = znew
  tmp[znew > 0] = NA
  lines(tmp,col=bcol)

  tmp = znew
  tmp[znew < 0] = NA
  lines(tmp,col=acol)
  
}
  



