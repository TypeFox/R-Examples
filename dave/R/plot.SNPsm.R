plot.SNPsm <-
function(x,...,out.seq=1,col=FALSE) {
   o.SNPsm<- x
   imax<- o.SNPsm$imax
   jmax<- o.SNPsm$jmax
   x<- o.SNPsm$sim.data
   tmap<- o.SNPsm$tmap
   igmap<- o.SNPsm$igmap # map of initial conditions
   cgmap<- matrix(as.character(igmap),nrow=imax)
   nt<- o.SNPsm$n.time.steps
   vegdef<- o.SNPsm$vegdef
   tsl<- o.SNPsm$time.step.length      # time step length
   frame<- o.SNPsm$frame
   out.seq<- max(out.seq,tsl)          # output sequence not shorter than tsl
   p.out<- as.integer(out.seq/tsl)     # output interval in time step
#
#  print map of initial conditions, igmap
#
   vegtypes<- c("1 Pinus","2 Carex","3 Festuca","4 Deschampsia","5 Trisetum","6 Aconitum")
   if(col == FALSE){
      par(mfrow=c(1,1),mar=c(0,0,0,0),oma=c(2,2,2,2))
      colors<- c(gray(0.85),gray(0.2),gray(0.45),gray(0.6),gray(0.7),gray(0.8))
      plot(c(1,jmax),c(1,imax),type="n",xlab="",ylab="",asp=1,axes=FALSE)
      for (i in 1:jmax) for (j in 1:imax) {
         if (frame[j,i]==1) points(i,(imax+1-j),cex=1.6,col=rgb(0,0,0,50,maxColorValue=255),pch=22)
         text(i,(imax+1-j),cgmap[j,i],cex=0.8,col=colors[igmap[j,i]],font=2,bty="n")
         }
      points(21,25,pch=1,cex=2.0,col="black")
      legend("right",c(1,2,3,4,5,6),vegtypes,cex=0.6,col=colors,bty="n")
   }
   if(col == TRUE) {
      colors<- c("darkolivegreen4","lightgreen","gold","darkorange","red1","darkred")
      plot(c(1,jmax),c(1,imax),type="n",xlab="",ylab="",asp=1,axes=FALSE)
      for (i in 1:jmax) for (j in 1:imax) {
         if (frame[j,i]==1) points(i,(imax+1-j),cex=1.6,col=rgb(0,0,0,50,maxColorValue=255),pch=22)
         text(i,(imax+1-j),cgmap[j,i],cex=0.8,col=colors[igmap[j,i]],font=2,bty="n")
         }
      points(21,25,pch=1,cex=2.0,col="black")
      legend("right",c(1,2,3,4,5,6),vegtypes,cex=0.6,col=colors,bty="n")
   }

#
#  Plot discrete map of state at t=1, using correlation
#
   par(mfrow=c(4,6),mar=c(0,0,0,0),omi=c(0,0,0,0),lwd=0.4)
   colors<- c(gray(0.10),gray(0.25),gray(0.40),gray(0.55),gray(0.70),gray(0.85))
   if(col == TRUE) colors<- c("darkolivegreen4","lightgreen","gold","darkorange","red1","darkred")
   for (i in 1:imax) for (j in 1:jmax) {
       tmap[i,j]<- which.max(cor(vegdef,x[i,j,1,1:6]))
   }
   plot(c(1,jmax),c(1,imax),type="n",xlab="",ylab="",asp=1,axes=FALSE)
   for (i in 1:jmax) for (j in 1:imax) {
      if (frame[j,i]!=1) points(i,(imax+1-j),pch=15,cex=0.8,col=colors[tmap[j,i]])
      if (frame[j,i]==1) points(i,(imax+1-j),pch=22,cex=0.8,col=rgb(0,0,0,100,maxColorValue=255))
   }
  legend("top","t= 0",xjust=0.5,box.lty=0,text.col="black",inset=c(0,-0.04))

   for (t in 1:nt) {
        if(t%%p.out == 0) {
#          cat("Output time =",t*tsl,"\n")
          for (i in 1:imax) for (j in 1:jmax) {
              tmap[i,j]<- which.max(cor(vegdef,x[i,j,t,1:6]))
          }
          plot(c(1,jmax),c(1,imax),type="n",xlab="",ylab="",asp=1,axes=FALSE)
          for (i in 1:jmax) for (j in 1:imax) {{
               if (frame[j,i]!=1) points(i,(imax+1-j),pch=15,cex=0.8,col=colors[tmap[j,i]])
               if (frame[j,i]==1) points(i,(imax+1-j),pch=22,cex=0.8,col=rgb(0,0,0,100,maxColorValue=255))
          }}
          legend("top",paste("t=",t*tsl),xjust=0.5,box.lty=0,text.col="black",inset=c(0,-0.04))

        } 
   }
  }
