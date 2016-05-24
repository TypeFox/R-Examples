plotCGP <-
function(object){
  plot_lim=c(min(object$Yp_jackknife,object$yobs)-0.01*(max(object$yobs)-min(object$yobs)),max(object$Yp_jackknife,object$yobs)+0.01*(max(object$yobs)-min(object$yobs)))
  plot(object$Yp_jackknife,object$yobs,type="p",pch=20,cex=0.8,xlim=plot_lim,ylim=plot_lim,xlab="Y Jackknife Predicted",ylab="Y Observed",main="Actually by Predicted Plot")
    
}
