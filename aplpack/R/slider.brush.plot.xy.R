slider.brush.plot.xy<-function(x,y=NULL,z=NULL,...)
{
  x.name<-deparse(substitute(x))
  y.name<-deparse(substitute(y))
  z.name<-deparse(substitute(z))
  if(length(x)<2) return("Error: x is of length 0 or 1")
  if(!is.null(y)){ 
    if(length(y)<2) return("Error: y must be a vector")
    if(length(x)!=length(y)) 
      return("Error: x and y must have the same length")
    x<-cbind(x,y) 
  }
  if(!is.null(z)){ 
    if(length(z)<2) return("Error: z must be a vector")
    if(length(x)!=length(z)) 
      return("Error: x and z must have the same length")
    x<-cbind(x,z) 
  }
  if(!is.matrix(x)&& !is.data.frame(x) && ncol(x)<3)
    return("Error: not enough variables")
  if("NULL"==y.name){x.name<-colnames(x)[1]; y.name<-colnames(x)[2]}
  if("NULL"==z.name){x.name<-colnames(x)[1]; z.name<-colnames(x)[3]}
  z<-x[,3]; y<-x[,2]; x<-x[,1]

  args<-list(...)
  if(!any("main"==names(args))) 
    args<-c(args,list(main=paste(x.name,"<-->",y.name)))
  if(!any("xlab"==names(args)))args<-c(args,list(xlab=x.name))
  if(!any("ylab"==names(args)))args<-c(args,list(ylab=y.name))
  do.call("plot.default",c(alist(x=x,y=y,pch=19),args))
  refresh<-function(...){
    zrange<-range(z); z1<-slider(no=1); z2<-slider(no=2)
    zmin<-z1; zmax<-z1+z2; ind<-zmin<=z&z<=zmax; pos<-par()$usr
    rect(pos[2],pos[4],pos[1]*.5+pos[2]*.5,pos[3]*.1+pos[4]*.9,
         col="white",border=NA)
    txt<-paste(z.name,"(red) in [",format(zmin,digits=4),",",
               format(zmax,digits=4),"]",sep="")
    text(pos[2],pos[4],txt,adj=c(1,1),col="red",cex=0.7)
    col<-c("black","red")[1+ind]
    points(x,y,col=col,pch=19,
           cex=if("cex" %in% names(args)) args$cex else 1)
  }
  z.min<-min(z); z.max<-max(z); delta<-(z.max-z.min)/100
  reset<-function(...){
    do.call("plot",c(alist(x=x,y=y,col="red",pch=19),args)); pos<-par()$usr #090216
    rect(pos[2],pos[4],pos[1]*.4+pos[2]*.6,pos[3]*.1+pos[4]*.9,
         col="white",border=NA)
  }
  slider(refresh,
         c("minimum of z","interval width of z"),
         c(z.min,0),c(z.max+delta,(z.max-z.min)+delta),
         c(delta,delta),c(z.min-delta,(z.max-z.min)/2),
         reset.function=reset
  )
  refresh()
  cat("use sliders to select interval for inking points\n")
}

