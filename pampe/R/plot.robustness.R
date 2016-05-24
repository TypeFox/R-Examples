plot.robustness <-
  function(x, ... ){
    
    robustness.object <- x
    
    
  if (class(robustness.object) != "robustness"){
    stop("Wrong object class")
  } 
  

  
  linewidth  <- matrix(1, 1, ncol(robustness.object)-2)
  linewidth <- append(linewidth, c(2,2), after = 0)
  matplot(robustness.object,
          type="l", xlab="", ylab="",
          col=c(1, 1, matrix("gray", 1, ncol(robustness.object)-2)),
          lty=c(1, 2, matrix(1, 1, ncol(robustness.object)-2)),
          lwd=linewidth, xaxt="n")
  ##Axis
  axis(1, at=c(seq(2, nrow(robustness.object), by=2)),
       labels=rownames(robustness.object)[c(seq(2, nrow(robustness.object), by=2))], las=3)
  
  title(main="Leave One Out Robustness Check")
  ##Legend
  legend("bottomleft",c("Actual", "Predicted", "leave-one-out"),
         col=c(1, 1, "gray"),lty=c(1,2,1),lwd=c(2,2, 1))
  
}

