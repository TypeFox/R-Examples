histplot <- 
  function(obj1,obj)
  {

# obj1 is the vector of the generated spikes
# obj is an "FPTdensity"-class object    

Ms <- length(obj1)    
N1t <- obj$time[length(obj$time)]
RStflag <- get("RStudioflag")

# produces histogram
if(!RStflag) dev.new()
hist(obj1,breaks=Ms/5,prob=TRUE,xlim=c(0,N1t),xlab = "time",main="Histogram of FP times")
par(new=TRUE)
plot(obj$time,obj$g0,type="l",lwd=1,axes=FALSE,main="", xlab = "", ylab = "",col='red')
legend("topright",lty=1, "FPT pdf by int",col="red")
}
