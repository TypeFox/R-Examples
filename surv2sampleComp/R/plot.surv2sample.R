######################################
# plot.surv2sample
######################################
plot.surv2sample=function(x, measure=NULL,baseline=0,...){

 if(is.null(measure)){
 	y=x$survfit ; class(y)="survfit" ; plot(y,...)
 }

 if(measure=="relative percentile"){
    xx=x$quanprobs
    if(baseline==0){y=x$contrast.ratio01}else{y=x$contrast.ratio10}
     yy=y[(nrow(y)-length(xx)+1):nrow(y),]
     plotCI(xx, yy[,1], uiw=yy[,3], liw=yy[,2], xlab="percent", ylab="relative time (95%CI)",...)
 }

}

