`msn.map` <-
function(msn,lat,...) {
  xa <- min(lat[,2])
  xb <- max(lat[,2])
  if ((xb-xa)<(360-(xb-xa))) {
    x1 <- max(xa-5,-180)
    x2 <- min(xb+5,180)
    y1 <- max(min(lat[,1])-5,-90)
    y2 <- min(max(lat[,1])+5,90)  
    map("world",xlim=c(x1,x2),ylim=c(y1,y2))
    points(lat[,2],lat[,1],...)
    for (i in 1:ncol(msn)) {
      for (j in 1:ncol(msn)) {
        if (msn[i,j]==1) {lines(c(lat[i,2],lat[j,2]),c(lat[i,1],lat[j,1]))}}}}
  else {
    nl <- lat[,2]
    for (i in 1:length(nl)) {
      if (nl[i]<0) nl[i]<-nl[i]+360}
    x1 <- max(min(nl)-5,0)
    x2 <- min(max(nl)+5,360)
    y1 <- max(min(lat[,1])-5,-90)
    y2 <- min(max(lat[,1])+5,90)  
    map("world2",xlim=c(x1,x2),ylim=c(y1,y2))
    points(nl,lat[,1],...)
    for (i in 1:ncol(msn)) {
      for (j in 1:ncol(msn)) {
        if (msn[i,j]==1) {lines(c(nl[i],nl[j]),c(lat[i,1],lat[j,1]))}}}}
}

