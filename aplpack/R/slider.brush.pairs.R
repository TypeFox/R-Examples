slider.brush.pairs<-function(x,...)
{
  args<-list(...)
  x.name<-deparse(substitute(x))
  if(missing(x)||length(x)<2) return("Error: x must be a vector")

  # preparation of data
  m<-dim(x)[2]; for(j in 1:m) x[,j]<-as.numeric(x[,j])
  mins<-apply(x,2,min); maxs<-apply(x,2,max)
  delta<-(maxs-mins)/100
  # initial plot
  varnames<-paste("var ",1:m,": ",colnames(x),sep="")
  dev.new(); par(mfrow=c(m,m),oma=c(0,0,0,0),mai=c(0,0,0,0),...)
  usr.array<-array(0,c(m,m,4)); axes<-FALSE
  for(i in 1:m){
    for(j in 1:m){
      # plot(x[,j],x[,i],axes=axes,type="p")
      do.call("plot",c(alist(x=x[,j],y=x[,i],type="p",axes=axes,xlab="",ylab=""),args))
      usr.array[i,j,] <- usr<-par()$usr
      if(i==j) text(usr[1],usr[4],varnames[i],adj=c(0,1),cex=5)
      rect(usr[1],usr[3],usr[2],usr[4])
    }
  }
  # update function
  refresh<-function(...){
    vmin<-slider(no=1)/100; vmax<-vmin+slider(no=2)/100
    vno <-slider(no=3)
    vmin<-mins[vno]*(1-vmin)+maxs[vno]*(vmin) 
    vmax<-mins[vno]*(1-vmax)+maxs[vno]*(vmax)
    ind <-vmin<=x[,vno] & x[,vno]<=vmax
    for(i in 1:m){
      for(j in 1:m){
        par(mfg=c(i,j),usr=usr.array[i,j,])
        points(x[ ,j],x[ ,i],col=0,cex=2,pch=19)
        points(x[ ind,j],x[ ind,i],col="red",pch=1)
        points(x[!ind,j],x[!ind,i],col="blue",pch=19)
      }
    }
  }
  # slider definition
  nt <- slider(refresh,
         c("lower limit (% of range)","width (% of range)", 
           paste("variable no: 1 ..",m)),
         c(0,0,1), c(100,100, m), c(1,1,1), c(0,30,1)
  )
  # tkwm.minsize(nt, "450", "150") # set width, height to prevent to small sizes
  refresh()
  cat("use sliders to select variable and interval width\n")
}

