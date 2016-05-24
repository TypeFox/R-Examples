screeArchetypal <- function(numArch,rss_lass_def,rss_step_ns,rss_step_alpha,rss_step_beta,ylim,main,
                            xlab,ylab,col=c("red","blue","green3"),axis2,seq,leg){
 archsPlot <- seq(length = numArch)
 plot(archsPlot, rss_lass_def, xaxt="n", yaxt="n", ylim=ylim, main=main, xlab=xlab, ylab=ylab, type="b") 
 points(archsPlot, rss_step_ns, type="b", col=col[1]) 
 points(archsPlot, rss_step_alpha, type="b", col=col[2], lty=2)
 points(archsPlot, rss_step_beta, type="b", col=col[3], lty=2)
 axis(1, at=archsPlot, labels=archsPlot)
 
 if(axis2){                          
  axis(2,at=seq,labels=seq)
 }
 if(leg){
  legend("topright",c("Archetypes",
                      expression(paste("Archetypoids from ", cand[ns])),
                      expression(paste("Archetypoids from ", cand[alpha])),
                      expression(paste("Archetypoids from ", cand[beta]))),
         lty = c(1,1,2,6),
         col=c("black",col[1],col[2],col[3]),
         text.col=c("black",col[1],col[2],col[3]))  
  }
}