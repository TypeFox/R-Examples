plot.penmodel <- function(x, ...){


 penplot(base.parms=x$parms.est[1:2], vbeta=x$parms.est[3:4], base.dist=attr(x,"base.dist") , variation="none", agemin=attr(x, "agemin"))

 data <- attr(x, "data")

 lines(survfit(Surv(time, status)~gender+mgene, data=data[data$proband==0 & !is.na(data$mgene),]), fun="event",     lty=c(2,1,2,1), col=c("red","red","blue","blue"))

  
}