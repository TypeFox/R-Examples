`complexity.ipec.CoxBoost` <-
function(response, x, boot.n.c=10, boost.steps=100, 
   eval.times=NULL, smooth=FALSE, full.data, ...){
   require(CoxBoost)
   require(locfit)
   require(survival)
   actual.data.c <- as.data.frame(x)
   xnames <- names(actual.data.c)
   time <- response[,"time"]
   status <- response[,"status"]
   actual.data.c$time <- time
   actual.data.c$status <- status
   
   n.c <- length(time)
   boot.index.c <- matrix(sapply(1:boot.n.c,function(b){sample(1:n.c,replace=TRUE)}), 
      byrow=TRUE, nrow=boot.n.c, ncol=n.c)
   not.in.sample.c <- list()
   for (i in 1:boot.n.c){ 
      not.in.sample.c[[i]] <- (1:n.c)[-unique(boot.index.c[i,])]
   }

   uncens <- which(status == 1)
   if (is.null(eval.times)){
      if(length(unique(time[uncens]))<100){
         w.eval.times.c <- c(0,sort(time[uncens]))
      } else {
      w.quantile <- 90
      space <- round(length(time[uncens])/100)
      index <- (1:w.quantile)*space
      w.eval.times.c <- c(0,sort(time[uncens])[index[index<length(time[uncens])]])
      }
   } else {
      w.eval.times.c <- sort(eval.times)
   }

   w.eval.times.c <- unique(w.eval.times.c)

   km.fit <- survfit(Surv(time, status)~1, data=actual.data.c)
   km.pred <- summary(object=km.fit, times=w.eval.times.c)$surv
   km.weight <- -1*diff(km.pred)

   fullcoxboost <- CoxBoost(time=time,status=status,x=x, stepno=boost.steps, ...)

   #class(fullcoxboost) <- c(class(fullcoxboost), "peperrinterinternal")

   w.full.apparent <- matrix(NA,nrow=boost.steps, ncol=length(w.eval.times.c))
   w.noinf.error <- matrix(NA,nrow=boost.steps, ncol=length(w.eval.times.c))

   for (m in 1:boost.steps){
      pe.w.full.apparent <- pmpec(object=fullcoxboost,
         data=actual.data.c, times=w.eval.times.c,
         model.args=list(complexity=m),
         external.time=full.data$time, external.status=full.data$status)
      w.full.apparent[m,] <- pe.w.full.apparent

      pe.w.noinf.error <- pmpec(object=fullcoxboost,
         data=actual.data.c, times=w.eval.times.c,
         model.args=list(complexity=m),
         external.time=full.data$time, external.status=full.data$status, type="NoInf")
      w.noinf.error[m,] <- pe.w.noinf.error
   }

   w.boot.error.wo <- array(dim=c(boost.steps, boot.n.c, length(w.eval.times.c)))

   for (actual.boot in 1:boot.n.c) {
      boot.fit <- CoxBoost(time=time[boot.index.c[actual.boot,]],
         status=status[boot.index.c[actual.boot,]], 
         x=x[boot.index.c[actual.boot,],], stepno=boost.steps, ...)

      #class(boot.fit) <- c(class(boot.fit), "peperrinterinternal")
    
      for (m in 1:boost.steps){
         w.pec.boot <- pmpec(object=boot.fit,
            data=actual.data.c[not.in.sample.c[[actual.boot]],], 
            times=w.eval.times.c, 
            model.args=list(complexity=m),
            external.time=full.data$time, external.status=full.data$status)
         w.boot.error.wo[m,actual.boot,] <- w.pec.boot
      }
   }

   w.mean.boot.error.wo <- apply(w.boot.error.wo,c(1,3),mean,na.rm=TRUE)
   w.boot632p.error.wo <- matrix(NA, nrow=boost.steps, ncol=length(w.eval.times.c))

   for (m in 1:boost.steps) {

      w.relative.overfit.wo <- ifelse(w.noinf.error[m,] > w.full.apparent[m,],
         (ifelse(w.mean.boot.error.wo[m,] < w.noinf.error[m,], 
         w.mean.boot.error.wo[m,], w.noinf.error[m,]) - w.full.apparent[m,])/
         (w.noinf.error[m,] - w.full.apparent[m,]), 0)
      w.weights <- .632/(1-.368*w.relative.overfit.wo)

      w.boot632p.error.wo[m,] <- (1-w.weights)*w.full.apparent[m,] + 
         w.weights*ifelse(w.mean.boot.error.wo[m,] < w.noinf.error[m,], 
         w.mean.boot.error.wo[m,], w.noinf.error[m,])
   }

   w.boot632p.error.smooth <- matrix(NA, ncol=length(w.eval.times.c), nrow=boost.steps)

   if (smooth==TRUE){
      for (i in 1:boost.steps){
         smoothdata <- as.data.frame(w.eval.times.c)
         smoothdata$error <- w.boot632p.error.wo[i,]
         smoother <- locfit(error~lp(w.eval.times.c), smoothdata)
         w.boot632p.error.smooth[i,] <- predict(object=smoother, newdata=w.eval.times.c)
      }
   } else {
      w.boot632p.error.smooth <- w.boot632p.error.wo
   }

   Lint.boot632p.error <- apply(t(w.boot632p.error.smooth[,1:(length(km.weight))])*km.weight,2,sum)

   min.ipec <- which.min(Lint.boot632p.error)
   min.ipec
}

