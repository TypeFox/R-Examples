plot.pampe <-
function(x, ... ){
  
  pampe.object <- x
  
  if (class(pampe.object) != "pampe"){
    stop("Wrong object class")
  } 
  
  if (length(pampe.object) > 3){
    par(ask=T)
  } 
   
  #Plot 1: Actual and counterfactual
  matplot(cbind(pampe.object$counterfactual[,1],pampe.object$counterfactual[,2]),
          type="l", ylab="", xlab="",
          col=1, lwd=2, xaxt="n")
  ##Axis labels & titles
  axis(1, at=c(seq(2, nrow(pampe.object$counterfactual), by=2)),
       labels=rownames(pampe.object$counterfactual)[c(seq(2, nrow(pampe.object$counterfactual), by=2))], las=3)
  ##Legend
  legend("bottomleft",c("Actual", "Predicted"),
         col=1, lty=c(1,2), lwd=2)
  ##Add a vertical line when the tr starts
  abline(v=length(pampe.object$model$fitted.values),lty=3, lwd=2)
  title(main="Actual and Counterfactual Path")
  
  #Plot 2: Control placebos
  if ("placebo.ctrl" %in% names(pampe.object)){
    par(ask=T)
    
    linewidth  <- matrix(2, 1, ncol(pampe.object$placebo.ctrl$mspe)-1)
    linewidth <- append(linewidth, 5, after = 0)
    #par(mai=c(2.2,1.2,0.15,0.3), mgp=c(1.6,0.5,0))
    matplot(pampe.object$placebo.ctrl$tr.effect, type="l", xlab="", ylab="", col=c("red",matrix(1, 1, ncol(pampe.object$placebo.ctrl$mspe)-1)), lty=c(1,matrix(2, 1, ncol(pampe.object$placebo.ctrl$mspe)-1)), lwd=linewidth, xaxt="n")
    axis(1, at=c(seq(2, nrow(pampe.object$counterfactual), by=2)),
         labels=rownames(pampe.object$counterfactual)[c(seq(2, nrow(pampe.object$counterfactual), by=2))], las=3)
    
    legend("bottomleft",c("Treated", "Controls"),col=c("red", 1),lty=c(1,2),lwd=c(5,2))
    abline(h=0,lty=3, lwd=3)
    abline(v=length(pampe.object$model$fitted.values),lty=3, lwd=3)
    title(main="Placebo Study. Control Reassignment.")
  } 
  
  
  
  #Plot 3: Control time
  if ("placebo.time" %in% names(pampe.object)){
    par(ask=T)
    
    mspe <- pampe.object$placebo.time$mspe
    linewidth  <- matrix(2, 1, ncol(mspe)-1)
    linewidth <- append(linewidth, 5, after = 0)
    matplot(pampe.object$placebo.time$tr.effect,
            type="l", xlab="", ylab="",
            col=c("red",matrix(1, 1, ncol(mspe)-1)),
            lty=c(1,matrix(2, 1, ncol(mspe)-1)),
            lwd=linewidth, xaxt="n")
    ##Axis
    axis(1, at=c(seq(2, nrow(pampe.object$counterfactual), by=2)),
         labels=rownames(pampe.object$counterfactual)[c(seq(2, nrow(pampe.object$counterfactual),
                                                            by=2))], las=3)
    
    ##Legend
    legend("bottomleft",c("Treated", "Time Controls"),
           col=c("red", 1),lty=c(1,2),lwd=c(5,2))
    ##Horizontal line
    abline(h=0,lty=3, lwd=2)
    abline(v=length(pampe.object$model$fitted.values),lty=3, lwd=3)
    title(main="Placebo Study. Time Reassignment.")
    
    
    
    
  } 
  
  
  
  
  
}

