#' plot2= flood level plot #2
#' @export
#' @param  obj is output object from the binteg function
#' @param  prob is how the probabilities are interpreted - assumed AEP 
plot2=function(obj,prob){
    
    # note - the 4th series, which allows for boundary effect to lower bound is not plotted
    if (prob=="ARI"){
    lt=c(2,1,1,2);
    cl=c("red","black","black","blue")
    plot(NA,NA,xlim=range(obj$p.ari),ylim=range(obj$zout),log="x",xlab="ARI",ylab="Design Vairable")
    abline(v=ceiling(obj$p.ari),col="grey",lty=2)
    lines(obj$p.ari,obj$zout[,1],lty=2,col="red",lwd=2)
    lines(obj$p.ari,obj$zout[,2],lty=1,col="black",lwd=2)
    lines(obj$p.ari,obj$zout[,4],lty=2,col="blue",lwd=2)
    lines(obj$p.ari,obj$zout[,3],lty=1,col="black",lwd=2)
   legend("topleft", c("Complete dependence","Specified dependence", "Independence"), 
            col = c("red","black","blue"),
            text.col = "green4", 
            lty = c(2,1,2), lwd=2, 
            merge = TRUE, bg = 'gray90')}

    if (prob=="AEP"){
    lt=c(2,1,1,2);
    cl=c("red","black","black","blue")
    plot(NA,NA,xlim=c(max(obj$p.aep),min(obj$p.aep)),ylim=range(obj$zout),log="x",xlab="AEP",ylab="Design Vairable")
    #abline(v=ceiling(obj$p.aep),col="grey",lty=2)
    lines(obj$p.aep,obj$zout[,1],lty=2,col="red",lwd=2)
    lines(obj$p.aep,obj$zout[,2],lty=1,col="black",lwd=2)
    lines(obj$p.aep,obj$zout[,4],lty=2,col="blue",lwd=2)
    lines(obj$p.aep,obj$zout[,3],lty=1,col="black",lwd=2)
    legend("topleft", c("Complete dependence","Specified dependence", "Complete independence"), 
            col = c("red","black","blue"),
            text.col = "green4", 
            lty = c(2,1,2), lwd=2, 
            merge = TRUE, bg = 'gray90')}
}
 
 
