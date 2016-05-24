`helio.plot` <-
function(c,cv=1,xvlab=c$xlab,yvlab=c$ylab,x.name="X Variables",y.name="Y Variables",lab.cex=1,wid.fact=0.75,main="Helio Plot",sub=paste("Canonical Variate",cv,sep=""),zero.rad=30,range.rad=20,name.padding=5,name.cex=1.5,axis.circ=c(-1,1),x.group=rep(0,dim(c$xstructcorr)[1]),y.group=rep(0,dim(c$ystructcorr)[1]),type="correlation"){
   #First, open up a new window
   plot.new()
   plot.window(c(-100,100),c(-100,100))   #Assume, for convenience, a +-100 world
   #Set the appropriate data, depending on whether this is a correlation or a variance plot
   if(type=="correlation"){
      xdat<-c$xstructcorr
      ydat<-c$ystructcorr
   }else if(type=="variance"){
      xdat<-c$xstructcorrsq
      ydat<-c$ystructcorrsq
   }else
      stop(paste("Plot type ",type," not supported.\n",sep=""))
   #Set radii for inner, middle, outer, and name circles
   ir<-zero.rad-range.rad
   mr<-zero.rad
   or<-zero.rad+range.rad
   nr<-zero.rad+range.rad+name.padding
   #Next, put the dividing line and any axis circles in place
   lines(c(0,0),c(-90,90))
   lines(mr*sin(2*pi*((0:100)/100)),mr*cos(2*pi*((0:100)/100)),lty=1)
   if(!is.null(axis.circ))
      for(i in 1:length(axis.circ))
         lines((mr+range.rad*axis.circ[i])*sin(2*pi*((0:100)/100)),(mr+range.rad*axis.circ[i])*cos(2*pi*((0:100)/100)),lty=3)
   #Label the two halves of the circle
   text(-50,95,label=x.name,cex=name.cex)
   text(50,95,label=y.name,cex=name.cex)
   #Label the ranges
   #text(rep(0,6),c(-45,-25,-5,5,25,45),label=c(1,0,-1,-1,0,1),pos=c(rep(2,3),rep(4,3)),cex=0.85,offset=0.1)
   #Get the number of variables in each set
   nx<-dim(xdat)[1]
   ny<-dim(ydat)[1]
   #Place rectangles and names for the x (left) variables
   for(i in 1:nx){
      #First, place rectangles
      if(xdat[i,cv]>0)   #Set fill color - black if positive, unfilled if negative
         bcol<-1
      else
         bcol<-NA
      bang<-(-pi/(nx+1))*i   #Determine the angle
      binc<-pi/(max(nx,ny)+1)*wid.fact/2
      bwinc<-ir*sin(binc)  #Determine the box width increment (based on binc at inner circle)
      bx<-vector()
      bx[1]<-mr*sin(bang)-bwinc*cos(-bang)
      bx[2]<-(mr+range.rad*xdat[i,cv])*sin(bang)-bwinc*cos(-bang)
      bx[3]<-(mr+range.rad*xdat[i,cv])*sin(bang)+bwinc*cos(-bang)
      bx[4]<-mr*sin(bang)+bwinc*cos(-bang)
      by<-vector()
      by[1]<-mr*cos(bang)-bwinc*sin(-bang)
      by[2]<-(mr+range.rad*xdat[i,cv])*cos(bang)-bwinc*sin(-bang)
      by[3]<-(mr+range.rad*xdat[i,cv])*cos(bang)+bwinc*sin(-bang)
      by[4]<-mr*cos(bang)+bwinc*sin(-bang)
      polygon(bx,by,col=bcol,lty=1)  #Draw the box
      #Next, place names
      text(nr*sin(bang),nr*cos(bang),label=xvlab[i],srt=(3*pi/2-bang)*(360/(2*pi)),pos=2,cex=lab.cex)
   }
   #Place rectangles and names for the y (right) variables
   for(i in 1:ny){
      #First, place rectangles
      if(ydat[i,cv]>0)   #Set fill color - black if positive, unfilled if negative
         bcol<-1
      else
         bcol<-NA
      bang<-(pi/(ny+1))*i   #Determine the angle
      binc<-pi/(max(nx,ny)+1)*wid.fact/2
      bwinc<-ir*sin(binc)  #Determine the box width increment (based on binc at inner circle)
      bx<-vector()
      bx[1]<-mr*sin(bang)-bwinc*cos(-bang)
      bx[2]<-(mr+range.rad*ydat[i,cv])*sin(bang)-bwinc*cos(-bang)
      bx[3]<-(mr+range.rad*ydat[i,cv])*sin(bang)+bwinc*cos(-bang)
      bx[4]<-mr*sin(bang)+bwinc*cos(-bang)
      by<-vector()
      by[1]<-mr*cos(bang)-bwinc*sin(-bang)
      by[2]<-(mr+range.rad*ydat[i,cv])*cos(bang)-bwinc*sin(-bang)
      by[3]<-(mr+range.rad*ydat[i,cv])*cos(bang)+bwinc*sin(-bang)
      by[4]<-mr*cos(bang)+bwinc*sin(-bang)
      polygon(bx,by,col=bcol,lty=1)  #Draw the box
      #Next, place names
      text(nr*sin(bang),nr*cos(bang),label=yvlab[i],srt=(pi/2-bang)*(360/(2*pi)),pos=4,cex=lab.cex)
   }
   #Perform grouping for the X variables, if needed.  0 means ungrouped, numbers above it are grouped.
   if((!is.null(x.group))&(max(x.group)>0)){
      for(i in unique(x.group))
         if(i>0){
            #Find first and last occurrence (they'd damn well better be sorted!)
            gvect<-(x.group%in%i)*(1:length(x.group))
            gvect<-gvect[gvect>0]
            minang<-min(gvect)*(-pi/(nx+1))
            maxang<-max(gvect)*(-pi/(nx+1))
            #Add a nice looking grouping thingee
             lines(((or+nr)/2)*sin((((0:100)/100)*(maxang-minang)+minang)),((or+nr)/2)*cos((((0:100)/100)*(maxang-minang)+minang)),lty=1)
            lines(c(((or+nr)/2)*sin(minang),nr*sin(minang)),c(((or+nr)/2)*cos(minang),nr*cos(minang)),lty=1)
            lines(c(((or+nr)/2)*sin(maxang),nr*sin(maxang)),c(((or+nr)/2)*cos(maxang),nr*cos(maxang)),lty=1)
         }
   }
   #Perform grouping for the Y variables, if needed.  0 means ungrouped, numbers above it are grouped.
   if((!is.null(y.group))&(max(y.group)>0)){
      for(i in unique(y.group))
         if(i>0){
            #Find first and last occurrence (they'd damn well better be sorted!)
            gvect<-(y.group%in%i)*(1:length(y.group))
            gvect<-gvect[gvect>0]
            minang<-min(gvect)*(pi/(ny+1))
            maxang<-max(gvect)*(pi/(ny+1))
            #Add a nice looking grouping thingee
             lines(((or+nr)/2)*sin((((0:100)/100)*(maxang-minang)+minang)),((or+nr)/2)*cos((((0:100)/100)*(maxang-minang)+minang)),lty=1)
            lines(c(((or+nr)/2)*sin(minang),nr*sin(minang)),c(((or+nr)/2)*cos(minang),nr*cos(minang)),lty=1)
            lines(c(((or+nr)/2)*sin(maxang),nr*sin(maxang)),c(((or+nr)/2)*cos(maxang),nr*cos(maxang)),lty=1)
         }
   }
   #Add title, if one is listed
   title(main=main,sub=sub)
}

