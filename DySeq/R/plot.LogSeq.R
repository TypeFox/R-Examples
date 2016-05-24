#'plot.LogSeq
#'
#'Generates interaction diagram for LogSeq Objects, see: \code{\link{LogSeq}}
#'
#'
#'@param x a LogSeq object, that should be printed.
#'@param y further arguments passed to or from other methods.
#'@param ... further arguments passed to or from other methods.
#'
#'@export

plot.LogSeq<-function(x, y, ...){
  
  if(length(x[[1]])==4){
    stop("Only one sequence was found!")
  }
             
    
  if(class(x)!="LogSeq") warning("x should be a LogSeq object!")

  lambdas<-x[[1]]

  intercept<-stats::t.test(lambdas[,1])
  actor<-stats::t.test(lambdas[,3])
  partner<-stats::t.test(lambdas[,2])
  interac<-stats::t.test(lambdas[,4])

  b<-c(intercept$estimate, actor$estimate, partner$estimate,interac$estimate)


  nn <- b[1]+(-1)*b[2]+(-1)*b[3]+b[4] # no prev. reaction
  yn <- b[1]+(1)*b[2]+(-1)*b[3]+(-1)*b[4] # only partner
  ny <- b[1]+(-1)*b[2]+(1)*b[3]+(-1)*b[4] # only actor
  yy <- b[1]+(1)*b[2]+(1)*b[3]+b[4] # both

  if(TRUE){   # yaxis=="prob" | yaxis="exp" for adding options later
    nn<-exp(nn)
    yn<-exp(yn)
    ny<-exp(ny)
    yy<-exp(yy)
  }

  if(TRUE){     # for yaxis=="prob"
    nn<-nn/(nn+1)
    yn<-yn/(yn+1)
    ny<-ny/(ny+1)
    yy<-yy/(yy+1)
  }

  plot(c(0, 1),c(yy, ny),
       xlim=c(-0.2,1.2), ylim=c(0,1),
       type="l",
       xlab="Partner reaction at t-1",
       ylab="Probability of actor reaction at timeinterval t",
       xaxt="n")
  axis(side = 1,
       at = c(0,1),
       labels = c("yes", "no"),
       tck=-.05)
  lines(c(0, 1),c(yn, nn), lty=2 )
  legend(0.5, .95, c("Actor reaction at t-1", "no actor reaction at t-1"), lty=1:2)



}


