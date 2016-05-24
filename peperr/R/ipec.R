ipec <- function(pe, eval.times, type=c("Riemann", "Lebesgue", "relativeLebesgue"), response=NULL){
   type <- match.arg(type)
   if (is.vector(pe)){
      pe <- matrix(pe, ncol=length(pe))
   }
   if (type == "Riemann"){
      ipec <- apply(t(pe[,1:(ncol(pe)-1), drop=FALSE])*diff(eval.times),2,sum)
   } else { 
      time <- response[,"time"]
      status <- response[,"status"]
      actual.data <- as.data.frame(matrix(1, ncol=1, nrow=length(time)))
      actual.data$time <- time
      actual.data$status <- status
      km.fit <- survfit(Surv(time, status)~1, data=actual.data)
      km.pred <- summary(object=km.fit, times=eval.times)$surv
      km.weight <- -1*diff(km.pred)
      if (type == "Lebesgue"){
         ipec <- apply(t(pe[,1:(length(km.weight)),drop=FALSE])*km.weight,2,sum)
      }
      if (type == "relativeLebesgue"){
         nullmodel <- aggregation.pmpec(full.data=actual.data, type="apparent", 
         response=response, x=matrix(1, ncol=1, nrow=length(time)), model=km.fit)
         relpe <- nullmodel - pe
         ipec <- apply(t(relpe[,1:(length(km.weight)),drop=FALSE])*km.weight,2,sum)
      }
   }
ipec
}

