spread.labels<-function (x,y,labels=NULL,ony=NA,offsets=NA,between=FALSE,
 linecol=par("fg"),srt=0,...) {

 if(missing(x)) 
  stop("Usage: spread.labels(x,y,labels,...)")
 nx<-length(x)
 ny<-length(y)
 if(is.na(ony)) {
  xylim<-par("usr")
  ony<-diff(range(x))/(xylim[2]-xylim[1])<diff(range(y))/(xylim[4]-xylim[3])
 }
 if(between) {
  if(length(linecol)==1) linecol<-rep(linecol,2)
  nlabels<-length(labels)
  if(ony) {
   newy<-seq(y[1],y[ny],length=nlabels)
   # put the labels in the middle
   labelx<-rep(mean(x),nlabels)
   # do the left lines
   segments(x[1:nlabels],y[1:nlabels],labelx-strwidth(labels)/2,newy,
    col=linecol[1])
   # now the right lines
   segments(x[(nlabels+1):ny],y[(nlabels+1):ny],labelx+strwidth(labels)/2,newy,
    col=linecol[2])
   text(labelx,newy,labels,srt=srt,...)
  }
  else {
   newx<-seq(x[1],x[nx],length=nlabels)
   # put the labels in the middle
   labely<-rep(mean(y),nlabels)
   # do the upper lines
   segments(x[1:nlabels],y[1:nlabels],newx+strheight(labels)/2,labely,
    col=linecol[1])
   # now the right lines
   #segments(x[(nlabels+1):ny],y[(nlabels+1):ny],labelx+strwidth(labels)/2,newy,
    #col=linecol[2])
   text(newx,labely,labels,srt=srt,...)
  }
 }
 else {
  if(ony) {
   sort.index<-sort.list(y)
   x<-x[sort.index]
   y<-y[sort.index]
   newy<-seq(y[1],y[ny],length=length(labels))
   if(is.na(offsets)) {
    offset<-diff(par("usr")[1:2])/4
    offsets<-rep(c(offset,-offset),ny/2+1)[1:ny]
   }
   segments(x+offsets,newy,x,y)
   text(x+offsets,newy,labels[sort.index],srt=srt,pos=c(4,2),...)
  }
  else {
   sort.index<-sort.list(x)
   x<-x[sort.index]
   y<-y[sort.index]
   nx<-length(x)
   newx <- seq(x[1],x[nx],length=length(labels))
   if(missing(offsets)) {
    offset<-diff(par("usr")[3:4])/4
    offsets<-rep(c(offset,-offset),nx/2+1)[1:nx]
   }
   text(newx,y+offsets,labels[sort.index],srt=srt,...)
   seggap<-strheight("M")*2
   if(all(offsets > 0) || all(offsets < 0))
    seggaps<-rep(seggap*(all(offsets)<0)*-1,length.out=nx)
   else seggaps<-rep(c(seggap,-seggap),length.out=nx)
   segments(newx,y+offsets-seggaps,x,y)
  }
 }
}
