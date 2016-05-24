simpleBoxplot<-function (x, param, orderGrp=NULL, file="boxplot_groups.pdf"){
    
    ## select measurements from array data list
    data <- select.measurements(x)
    pdf(file=file)   
    for ( i in 1:ncol(data[[1]])){
      if(is.null(orderGrp)){
        arrayDat <- data.frame(expr=data[[1]][,i], param=data[[4]][,param])
      }else{
        arrayDat <- data.frame(expr=data[[1]][,i], param=factor(data[[4]][,param],levels=orderGrp))
      }
      boxplot(expr ~ param, data=arrayDat, main=c("target: ",data[[3]]["target",i]),las=2,ylab="signal intensity [a.u.]",xlab=param)
      stripchart(expr ~ param, data=arrayDat, add=T,vertical=T,method="jitter",jitter=0.3,col="red")
    }
    dev.off()
}
