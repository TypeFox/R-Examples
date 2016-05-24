 root.plot <- function(A){

  L <- eigen(A)$values
  nL <- length(L)

  # determine xlim and ylim to include the origin
  xmin <- min(Re(L)); xmax <- max(Re(L))
   if(xmin<0 && xmax<0) xmax <- 0
   if(xmin>0 && xmax>0) xmin <- 0
   if(xmin==0 && xmax==0) {xmin<- -1;xmax <- 1}
  ymin <- min(Im(L)); ymax <- max(Im(L))
   if(ymin<0 && ymax<0) ymax <- 0
   if(ymin>0 && ymax>0) ymin <- 0
   if(ymin==0 && ymax==0) {ymin<- -1;ymax <- 1}

   plot(Re(L[1]),Im(L[1]),xlab="Re",ylab="Im", xlim=c(xmin,xmax),ylim=c(ymin,ymax))
   for(i in 2:nL){
   if(c(Re(L[i]),Im(L[i])) - c(Re(L[i-1]),Im(L[i-1])) <= c(0,0)){
    dx <- 0.05*xmax; dy <- 0.05*ymax
   } else {dx <- 0; dy<-0}
     x <- Re(L[i])+dx
     y <- Im(L[i])+dy
   lines(x,y,type="p")
  }
  abline(h=0,v=0, col="grey")
} # end of function

