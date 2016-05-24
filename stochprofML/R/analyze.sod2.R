analyze.sod2 <-
function(model="LN-LN",TY=2,use.constraints=F) {
# Estimates the chosen model for the SOD2 dataset provided in this package.
#
# - model can be "LN-LN", "rLN-LN" or "EXP-LN"
# - TY is the number of types of cells that is assumed in the model
   
   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)
   d.sum.of.mixtures <- NULL
   sod2 <- NULL  
   rm(d.sum.of.mixtures)   
   rm(sod2)
   
   # model has to be defined
   if (!(model %in% c("LN-LN","rLN-LN","EXP-LN"))) {
      stop("Unknown model.")
   }   
   
   # settings
   dataset <- as.matrix(sod2,ncol=1)
   m <- 1
   
   # estimate
   result <- stochprof.loop(model=model,dataset=dataset,n=10,TY=TY,genenames="SOD2",fix.mu=F,loops=10,until.convergence=F,print.output=F,show.plots=T,plot.title="SOD2",use.constraints=use.constraints)
   par(ask=T)
   
   
   # histogram of data and pdf of estimated model
   hist(dataset,main="SOD2",breaks=20,xlab="Sum of mixtures of lognormals",ylab="Density",freq=F,col="lightgrey")
   x <- seq(round(min(dataset)),round(max(dataset)),(round(max(dataset))-round(min(dataset)))/100)
   # theoretical density
   mle <- result$mle
   # estimated p
   if (TY>1) {
      p.est <- mle[1:(TY-1)]
      p.est <- c(p.est,1-sum(p.est))
   }
   else {
      p.est <- 1
   }
   # estimated mu, sigma and lambda (if applicable);
   # calculate density
   if (model=="LN-LN") {
      mu.est <- mle[TY:(length(mle)-1)]
      sigma.est <- rep(mle[length(mle)],TY)
      
      y <- d.sum.of.mixtures(x,n=10,p.vector=p.est,mu.vector=mu.est,sigma.vector=sigma.est,logdens=F)
   }
   else if (model=="rLN-LN") {
      mu.est <- mle[TY:((m+1)*TY-1)]
      sigma.est <- mle[((m+1)*TY):length(mle)]
      
      y <- d.sum.of.mixtures(x,n=10,p.vector=p.est,mu.vector=mu.est,sigma.vector=sigma.est,logdens=F)
   }
   else if (model=="EXP-LN") {
      if (TY>1) {
         mu.est <- mle[TY:(2*(TY-1))]
         sigma.est <- rep(mle[2*(TY-1)+1],TY-1)
      }
      else {
         mu.est <- NULL
         sigma.est <- NULL
      }
      lambda.est <- mle[length(mle)]   
      
      y <- d.sum.of.mixtures(x,n=10,p.vector=p.est,mu.vector=mu.est,sigma.vector=sigma.est,lambda=lambda.est,logdens=F)
   }
   
   lines(x,y,col="blue",lwd=3)
   par(ask=F)   
   
   return(invisible(result))
}
