basemap <-
function(xlim,ylim,xlab='Longitude',ylab='Latitude',
                    bg='lightblue',...){
  # more-or-less appropriateaspect ratio
  asp <- 1/cos(sum(ylim)*pi/360) 
  plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,asp=asp,...)
  # background colour (does not work properly if you re-size plot window)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = bg)
}

