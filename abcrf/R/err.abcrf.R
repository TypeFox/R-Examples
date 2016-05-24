err.abcrf <- function(x)
{
  rf <- x$model.rf
  modindex <- rf$y
  nmod <- nlevels(modindex)
  sumsta <- x$sumsta
  mimi <- predict(rf,sumsta,predict.all=TRUE)$individual
  if (rf$ntree<10) stop("the number of trees in the forest should be greater than 10")
  sequo <- seq(10,rf$ntree,by=10)
  if (is.null(rf$err.rate))
  {
    res <- rep(0,length(sequo))
    h <- 0
    pb=txtProgressBar(min=0,max=length(sequo),style=3)
    for (j in sequo) 
    {
      mama <- sample(paste(1:nmod),nrow(sumsta),replace=TRUE)
      for (i in 1:nrow(sumsta)) 
      {
        outbag <- (1:j)[rf$inbag[i,1:j]==0]
        if (length(outbag)>0) mama[i] <- names(which.max(table(as.factor(mimi[i,outbag]))))
      }
    h <- h+1
    res[h] <- mean(as.factor(mama)!=modindex)
    setTxtProgressBar(pb,h)
    }
    close(pb)
  } else res <- rf$err.rate[sequo,1]
  plot(sequo,res,ylab="Prior error rate",xlab="Number of trees",type="l")
  cbind(ntree=sequo,error.rate=res)
}
