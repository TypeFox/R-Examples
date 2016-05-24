slider.density<-function(x,panel=rug,...)
{
  x.name<-deparse(substitute(x))
  if(missing(x)||length(x)<2) return("Error: x must be a vector")

  args<-list(...)
  if(!any("main"==names(args))) args<-c(args,list(main=x.name))
  kernel<-c("gaussian",  "epanechnikov","rectangular",
            "triangular","biweight",    "cosine",     "optcosine")
  slider(obj.name="kno",obj.value=1)
  refresh<-function(...){
    width<-slider(no=1)*diff(range(x))/100
    kno<-slider(obj.name="kno"); kernel<-kernel[kno]
    xy<-density(x,width=width,kernel=kernel)
    do.call("plot",c(alist(x=xy),args))
    title(paste("\n\nwidth =",signif(width,4),", kernel =",kernel))
    panel(x)
  }
  set.kernel<-function(...){  
    kernel<-slider(no=2)
    slider(obj.name="kno",obj.value=kernel)
    refresh()
  }
  bw.default<-diff(range(x))/density(x)$bw
  nt <- slider(c(refresh,set.kernel),
         c("width (% of range)","kernel"),
         c(.1,1),c(100,7),c(.1,1),c(bw.default,1)
  )
  # tkwm.minsize(nt, "300", "110") # set width, height to prevent to small sizes
  refresh()
  cat("use slider to select width of window and to select kernel:\n")
  print(cbind("no"=1:7,"kernel"=kernel))
}

