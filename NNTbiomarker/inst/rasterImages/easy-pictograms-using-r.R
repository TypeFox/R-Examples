# http://www.r-bloggers.com/easy-pictograms-using-r/
# requires image to be read in by readPNG or similar and supplied as "icon"
# To do: allow for non-integer n
pictogram<-function(icon,n,grouplabels="",
                    hicons=20,vspace=0.5,labprop=0.2,labelcex=1) {
  if(is.list(icon)) {
    licon<-icon
  } else {
    licon<-list(icon)
    for (i in 2:length(n)) {
      licon[[i]]<-icon
    }
  }
  library(reshape)
  sumn<-sum(n)
  group<-untable(df=matrix((1:length(n)),ncol=1),num=n)
  vicons<-ceiling(n/hicons)
  allv<-sum(vicons)
  tail<-n%%hicons
  # dim[1] is the height, dim[2] the width:
  devaspect<-dev.size(units="px")[1]/dev.size(units="px")[2]
  xlength<-1
  # get dims of all elements of licon, find greatest aspect and set ylength
  getdim<-function(z) {
    aspect<-dim(z)[1]/dim(z)[2]
    return(aspect)
  }
  all.ylengths<-unlist(lapply(licon,getdim))
  ylength<-max(all.ylengths)
  all.ylengths<-untable(df=matrix(all.ylengths,ncol=1),num=n)
  ytop<-allv*ylength
  if(devaspect*hicons<allv) warning("Icons may extend above the top of the graph")
  # vector of icons per row
  iconrow<-as.vector(as.matrix(rbind(rep(hicons,length(vicons)),tail)))
  # vector for how many times to repeat elements of iconrow
  reprow<-as.vector(as.matrix(rbind((vicons-1),rep(1,length(vicons)))))
  perrow<-untable(df=matrix(iconrow,ncol=1),num=reprow)
  spacing<-NULL
  for (i in 1:(length(n))) {
    spacing<-c(spacing,rep((i-1)*vspace*ylength,n[i]))
  }
  y0<-spacing+(ylength*untable(df=matrix((1:allv)-1,ncol=1),num=perrow))
  y1<-y0+all.ylengths
  # there are more elegant ways to make x0, but for now...
  x0<-NULL
  for (i in 1:(length(perrow))) {
    x0<-c(x0,(0:(perrow[i]-1)))
  }
  x1<-x0+xlength
  leftplot<-floor(-(labprop*hicons))
  plot(c(leftplot,hicons),c(0,(devaspect*hicons)),
       type="n",bty="n",ylab="",xlab="",xaxt="n",yaxt="n")
  lines(x=c(0,0),y=c(min(y0)-(ylength/2),max(y1)+(ylength/2)))
  for (i in 1:sumn) {
    rasterImage(image=licon[[group[i]]],xleft=x0[i],xright=x1[i],
                ytop=y1[i],ybottom=y0[i])
  }
  # find positions for labels
  ylabpos<-rep(NA,length(n))
  for (i in 1:length(n)) {
    ylabpos[i]<-(max(y1[group==i])+min(y0[group==i]))/2
  }
  text(x=leftplot/2,y=ylabpos,labels=grouplabels,cex=labelcex)
}
