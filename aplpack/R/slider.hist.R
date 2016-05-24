slider.hist<-function(x,panel=rug,...)
{
  x.name<-deparse(substitute(x))
  if(missing(x)||length(x)<2) return("Error: x must be a vector")

  args<-list(...)
  args<-args[names(args)!="breaks"]
  ClassNumber<-length(hist(x,plot=FALSE)$breaks)
  if(!any("main"==names(args)))args<-c(args,list(main=x.name))
  refresh<-function(...){
    xrange<-range(x); num<-slider(no=1)
    breaks<-seq(xrange[1],xrange[2],length=num+1)
    do.call("hist",c(alist(x=x,breaks=breaks),args))
    panel(x)
  }
  slider(refresh,"ClassNumber",1,100,1,ClassNumber);  refresh()
  "use slider to select number of classes"
}

