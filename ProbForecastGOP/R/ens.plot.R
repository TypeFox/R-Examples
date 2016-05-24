"ens.plot" <-
function(grid,lims,x.lim,y.lim,title){  
   size=65  
   hor.crd <- seq(x.lim[1],x.lim[2],length = size)  
   ver.crd <- seq(y.lim[1],y.lim[2],length = size)  
   n.level <- 5
   US.map <- 0
   image.plot(hor.crd, ver.crd, grid, xlim=c(min(hor.crd), max(hor.crd)), ylim = c(min(ver.crd),max(ver.crd)),zlim=lims,legend.width=.015,legend.shrink=.8,
   main=title,xlab="",ylab="",offset = 0.02,col=rainbow(100,start=0,end=0.85)[100:1])
   if(min(hor.crd) <= -124.7 & max(hor.crd) >= -124.7){
    US.map <- 1
   }
   if(min(hor.crd) >= -124.7 & max(hor.crd) <= -67.1){
    US.map <- 1
   }
   if(min(hor.crd) <= -67.1 & max(hor.crd) >= -67.1){
    US.map <- 1
   }
   if(min(ver.crd) <= 25.2 & max(ver.crd) >= 25.2){
    US.map <- 1 + US.map
   }
   if(min(ver.crd) >= 25.2 & max(ver.crd) <= 49.4){
    US.map <- 1 + US.map
   }
   if(min(ver.crd) <= 49.4 & max(ver.crd) >= 49.4){
    US.map <- 1 + US.map
   }
   if(US.map==2){
     US(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
     world(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)     
   }  
   
   if(US.map < 2){
       world(xlim=c(min(hor.crd), max(hor.crd)), ylim =c(min(ver.crd), max(ver.crd)),add=TRUE,col=1,lwd=2)
   }
}
