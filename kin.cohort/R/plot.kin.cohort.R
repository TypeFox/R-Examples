`plot.kin.cohort` <-
function(x,what=c("cr", "hr", "crr"),min.age=min(x$knots), max.age=max(x$knots), max.y, type, add=FALSE, color, line, ...){
   if(missing(color)) color<-rep(0,x$ngeno.rel)
   else if(length(color)< x$ngeno.rel) color<-rep(color,x$ngeno.rel)
   if(missing(line)) line<-rep(0,x$ngeno.rel)
   else if(length(line)< x$ngeno.rel) line<-rep(line,x$ngeno.rel)

   if (what[1]=="cr"){
      if (missing(max.y)) max.y=max(x$cumrisk[,1:x$ngeno.rel])
      if (missing(type)) type="l"
      min.risk<-pmax(0,x$cumrisk[match(min.age,x$knots),],na.rm=TRUE)

      if (add) 
      lines (c(min.age,x$knots),c(min.risk[1],x$cumrisk[,1]),col=1+color[1],type=type,lwd=2+line[1], ...)
      else
      plot (c(min.age,x$knots),c(min.risk[1],x$cumrisk[,1]),col=1+color[1],type=type,xlim=c(0,max.age),ylim=c(0,max.y),xlab="Age",ylab="Cumulative Risk",lwd=2+line[1], ...)
      for (i in 2:x$ngeno.rel)
         lines(c(min.age,x$knots),c(min.risk[i],x$cumrisk[,i]),col=i+color[i],type=type,lwd=1+line[i], ...)
      if (!add) 
      legend(0,max.y,colnames(x$cumrisk)[x$ngeno.rel:1],fill=(x$ngeno.rel:1+color))
   }
   if (what[1]=="hr" ){
      if(!inherits(x,"chatterjee")) stop("Hazard ratio not available")
      if (missing(type)) type="s"
      if (missing(max.y)) max.y=max(c(x$hazard[,(x$ngeno.rel+1)],1/x$hazard[,(x$ngeno.rel+1)]))
      if (add) 
      lines (c(min.age,x$knots),c(0,x$hazard[,(x$ngeno.rel+1)]),col=2+color[1],type=type,lwd=1+line[1], ...)
      else
      plot (c(min.age,x$knots),c(1,x$hazard[,(x$ngeno.rel+1)]),col=2+color[1],type=type,xlim=c(0,max.age),ylim=c(1/max.y,max.y),xlab="Age",ylab="Hazard Ratio",log="y",lwd=1+line[1], ...)
      for (i in 2:(x$ngeno.rel-1))
          if (i>2)
             lines (c(min.age,x$knots),c(1,x$hazard[,(x$ngeno.rel+i)]),col=(i+1)+color[i],type=type,lwd=1+line[i],...)
      lines(c(0,max(x$knots)),c(1,1),col=1+color,lty=1)
    }
   if (what[1]=="crr"){
      if (missing(type)) type="s"
      if (missing(max.y)) max.y=max(c(x$cumrisk[,(x$ngeno.rel+1)],1/x$cumrisk[,(x$ngeno.rel+1)]))
      if (add) 
      lines (c(min.age,x$knots),c(0,x$cumrisk[,(x$ngeno.rel+1)]),col=1+color[1],type=type,lwd=1+line[1], ...)
      else
      plot (c(min.age,x$knots),c(1,x$cumrisk[,(x$ngeno.rel+1)]),col=1+color[1],type=type,xlim=c(0,max.age),ylim=c(1/max.y,max.y),xlab="Age",ylab="Cumulative Risk Ratio",log="y",lwd=1+line[1], ...)
      for (i in 2:(x$ngeno.rel-1))
          if (i>2)
             lines (c(min.age,x$knots),c(1,x$cumrisk[,(x$ngeno.rel+i)]),col=i+color[i],type=type,lwd=1+line[i],...)
      abline(h=1,col=4+color,lty=3)
   }
}
