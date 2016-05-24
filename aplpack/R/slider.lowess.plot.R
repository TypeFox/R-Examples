slider.lowess.plot<-function(x,y=NULL,...)
{
  # slider function to draw lowess smoother, pwolf 080525
  x.name<-deparse(substitute(x))
  y.name<-deparse(substitute(y))
  if(length(x)<2) return("Error: x is of length 0 or 1")
  if(!is.null(y)){ 
    if(length(y)<2) return("Error: y must be a vector")
    if(length(x)!=length(y)) 
      return("Error: x and y must have the same length")
    x<-cbind(x,y) 
  }
  if(!is.matrix(x) && !is.data.frame(x)){ 
    x<-cbind(seq(x),x) 
    y.name<-x.name; x.name<-"index"
  }
  if(is.null(y.name)){x.name<-colnames(x)[1]; y.name<-colnames(x)[2]}
  y<-x[,2]; x<-x[,1]

  args<-list(...)
  refresh<-function(...){
    f<-slider(no=1)
    iter<-slider(no=2)
    xy<-lowess(x,y,f=f,iter=iter)
    # plot(x,y,bty="n")
    do.call("plot",c(alist(x,y,bty="n"),args))
    lines(xy)
    title(paste("\n\nlowess: f =",signif(f,4),", iterations =",iter))
    lines(xy)
  }
  slider(refresh,
         c("smoother span","iterations"),
         c(.01,1),c(1,7),c(.01,1),c(2/3,3)
  )
  refresh()
  cat("use slider to select smoother span!\n")
}

