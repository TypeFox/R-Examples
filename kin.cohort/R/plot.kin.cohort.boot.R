`plot.kin.cohort.boot` <-
function(x, conf = 0.95, what=c("cr","hr","crr"),min.age=min(x$knots), max.age=max(x$knots), age.start=0, start.ref, max.y, type, median=FALSE, add=FALSE, color, line, ...){
#
# only works for 1 gene
#
  if(missing(color)) color<-rep(0,2)
  else if(length(color)==1) color<-rep(color,2)
  if(missing(line)) line<-rep(0,3)
  else if(length(line)<3) line<-rep(line,3)
  
  if (!(what[1]%in%  c("cr","hr","crr"))) stop("valid args for what:  'cr' 'hr' 'crr'")
  if (!inherits(x,"print.kin.cohort.boot")) x<- print.kin.cohort.boot(x, conf=conf, show=FALSE)

  if( x$ncateg==4) stop("Sorry, only 1 gene can be plotted with bootstrap CI. This will be fixed ASAP", call.=FALSE)

  if (what[1]=="cr"){
      if (missing(max.y)) max.y<-max(x$cumrisk[,1:8])
      if (missing(type)) type<-"l"
      min.risk<-pmax(0,x$cumrisk[match(min.age,x$knots),],na.rm=TRUE)
#     lines (c(min.age,x$knots),c(min.risk[1],x$cumrisk[,1]),col=1+color[1],type=type, ...)

      if (add) 
      lines(c(min.age,x$knots),c(min.risk[1+median],x$cumrisk[,1+median]),col=1+color[1],type=type,lty=1,lwd=2+line[1], ...)
      else
      plot (c(min.age,x$knots),c(min.risk[1+median],x$cumrisk[,1+median]),col=1+color[1],type=type,lty=1,lwd=2+line[1], xlim=c(age.start,max.age),ylim=c(0,max.y),xlab="Age",ylab="Cumulative Risk", ...)
      lines(c(min.age,x$knots),c(min.risk[3],x$cumrisk[,3]),col=1+color[1],type=type,lty=3,lwd=1+line[3] )
      lines(c(min.age,x$knots),c(min.risk[4],x$cumrisk[,4]),col=1+color[1],type=type,lty=3,lwd=1+line[3] )
      lines(c(min.age,x$knots),c(min.risk[5+median],x$cumrisk[,5+median]),col=2+color[2],type=type,lty=1,lwd=2+line[2] )
      lines(c(min.age,x$knots),c(min.risk[7],x$cumrisk[,7]),col=2+color[2],type=type,lty=3,lwd=1+line[3] )
      lines(c(min.age,x$knots),c(min.risk[8],x$cumrisk[,8]),col=2+color[2],type=type,lty=3,lwd=1+line[3] )
      if(!add)
      legend(0,max.y,colnames(x$cumrisk)[c(5,1)],fill=c(2+color[2],1+color[1]))
   }
   if (what[1]=="hr"){
      if(!inherits(x,"chatterjee")) stop("Hazard Ratio not available")

      if (missing(type)) type<-"s"
      if (missing(max.y)) {
         val<-c(x$hazard[,9:12],1/x$hazard[,9:12])
         max.y<-max(val[is.finite(val)])
      }
      if(missing(start.ref))start.ref<-x$hazard[1,9+median]
      if(add)
      lines(c(min.age,x$knots),c(start.ref,x$hazard[,9+median]),type=type,col=1+color[1],lwd=1+line[1], ...)
      else
      plot (c(min.age,x$knots),c(start.ref,x$hazard[,9+median]),type=type,col=1+color[1],lwd=1+line[1], xlim=c(0,max.age),ylim=c(1/max.y,max.y),xlab="Age",ylab="Hazard Ratio",log="y", ...)
      lines(c(min.age,x$knots),c(0,x$hazard[,11]),col=2+color[2],type=type,lty=3,lwd=1+line[2] )
      lines(c(min.age,x$knots),c(0,x$hazard[,12]),col=2+color[2],type=type,lty=3,lwd=1+line[2] )
      if(!add)
      abline(h=1,col=4,lty=3)
   }
   if (what[1]=="crr"){

      if (missing(type)) type="s"
      if (missing(max.y)) {
         val<-c(x$cumrisk[,9:12],1/x$cumrisk[,9:12])
         max.y<-max(val[is.finite(val)])
      }
      if(missing(start.ref))start.ref<-x$cumrisk[1,9+median]
      if(add)
      lines(c(min.age,x$knots),c(start.ref,x$cumrisk[,9+median]),type=type,col=1+color[1],lwd=1+line[1], ...)
      else
      plot (c(min.age,x$knots),c(start.ref,x$cumrisk[,9+median]),type=type,col=1+color[1],lwd=1+line[1],xlim=c(0,max.age),ylim=c(1/max.y,max.y),xlab="Age",ylab="Cumulative Risk Ratio",log="y", ...)
      lines(c(min.age,x$knots),c(0,x$cumrisk[,11]),col=2+color[2],type=type,lty=3,lwd=1+line[2] )
      lines(c(min.age,x$knots),c(0,x$cumrisk[,12]),col=2+color[2],type=type,lty=3,lwd=1+line[2] )
      if(!add)
      abline(h=1,col=4,lty=3)
   }}
