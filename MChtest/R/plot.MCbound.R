"plot.MCbound" <-
function(x, rdigit=4, plimit=500, ...){
    op<-par(mfrow=c(2,1),oma=c(0,0,3,0))
   plot(c(0,max(x$N)),range(c(0,x$S)),type="n",axes=FALSE,
       main="Stopping boundary",xlab="N=number of replications",ylab="S=number of successes",...)
   if (length(x$p.value)>plimit){
       lines(x$N,x$S,lty=1)
   }
   else{ points(x$N,x$S,pch=20,cex=.5) }

    axis(1)
    axis(2)
   
    YLIM<- max(c(-x$ci.lower+x$p.value,x$ci.upper-x$p.value)) 
    plot(c(0,1),
      c(-YLIM,YLIM),
      xlab="estimated p value",ylab="ci.limits-estimated p value",type="n",
      main=paste(100*x$conf.level,"% Confidence Limits on P values"),...)
    if (length(x$p.value)>plimit){
    lines(x$p.value,x$p.value-x$p.value,lty=2)
    lines(x$p.value,x$ci.lower-x$p.value,lty=1)
    lines(x$p.value,x$ci.upper-x$p.value,lty=1)
    }
    else{
    points(x$p.value,x$p.value-x$p.value,pch=20)
    points(x$p.value,x$ci.lower-x$p.value,pch=21)
    points(x$p.value,x$ci.upper-x$p.value,pch=21)
    }

    #plot(c(0,1),
    #  c(0,1),
    #  xlab="p value",ylab="p value + or - confidence limits",type="n",
    #  main=paste("Plot of ",x$Name))
    #lines(x$p.value,x$p.value,lty=1)
    #lines(x$p.value,x$ci.lower,lty=2)
    #lines(x$p.value,x$ci.upper,lty=2)



    TITLE<-paste("Monte Carlo stopping boundary of type=",x$type,"\n with")
    for (i in 1:length(x$parms)){
        TITLE<-paste(TITLE,names(x$parms)[i],"=",round(x$parms[i],rdigit))
    }
    mtext(TITLE,outer=TRUE)
 
    par(op)
}

