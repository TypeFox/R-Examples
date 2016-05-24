setMethod("plot",signature=signature("OPT"),
    definition=function(x) {
         obj<-x@Par
         LB<-obj@LB
         UB<-obj@UB
         gr<-obj@grid
         x1<-seq(LB,UB,gr)
         ds<-obj@ds
         if(obj@fid=="MOPT") plot(x1,ds,cex=.3,main="Verify the Multi-obj optimal design",ylab="Sensitivity function",xlab="Log Dose")
         if(obj@fid=="ceff1") plot(x1,ds,cex=.1,main="Verify the c-optimal design for ED50",ylab="Sensitivity function",xlab="Log Dose")
         if(obj@fid=="ceff2") plot(x1,ds,cex=.3,main="Verify the c-optimal design for MED",ylab="Sensitivity function",xlab="Log Dose")
         if(obj@fid=="Deff") plot(x1,ds,cex=.3,main="Verify the D-optimal design",ylab="Sensitivity function",xlab="Log Dose")
    }
)

      