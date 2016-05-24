generate.toydata <-
function(model="LN-LN") {
# Generates a synthetic dataset (without measurement error) for one gene 
# from the chosen model and estimates the parameters from this data.
#
# model can be "LN-LN", "rLN-LN" or "EXP-LN".


   # definition of variables (necessary for CMD R check)
   # (these variables will be initialized later, but they are not visible as global functions/data)
   r.sum.of.mixtures <- NULL
   d.sum.of.mixtures <- NULL
   rm(r.sum.of.mixtures)
   rm(d.sum.of.mixtures)


   # model has to be defined
   if (!(model %in% c("LN-LN","rLN-LN","EXP-LN"))) {
      stop("Unknown model.")
   }   
   set.model.functions(model)

   ############
   # settings #
   ############
   
   # number of types
   TY <- 2
   # probabilities for both types
   p.vector <- c(0.2,0.8)
   # mean and standard deviation
   if (model=="LN-LN") {
      # log-mean
      mu.vector <- c(1.5,-1.5)
      # log-standard deviation
      sigma.vector <- rep(0.2,2)
   }
   else if (model=="rLN-LN") {
      # log-mean
      mu.vector <- c(1.5,-1.5)
      # log-standard deviation
      sigma.vector <- c(0.2,0.6)   
   }   
   else if (model=="EXP-LN") {
      # log-mean
      mu.vector <- c(1.5)
      # log-standard deviation
      sigma.vector <- rep(0.2,TY-1)
      # exponential rate
      lambda.vector <- 0.5   
   }   
   # number of independent samples (tissue samples)
   k <- 500
   # number of cells in each tissue sample
   n <- 10

   ##########################
   # generate and plot data #
   ##########################
   
   # generate
   m <- 1
   if (model %in% c("LN-LN","rLN-LN")) {
      dataset <- r.sum.of.mixtures(k=k,n=n,p.vector=p.vector,mu.vector=mu.vector,sigma.vector=sigma.vector)
      xlabel <- "Sum of mixtures of lognormals"
   }
   else if (model=="EXP-LN") {
      dataset <- r.sum.of.mixtures(k=k,n=n,p.vector=p.vector,mu.vector=mu.vector,sigma.vector=sigma.vector,lambda=lambda.vector) 
      xlabel <- "Sum of mixtures of lognormals and exponentials"  
   }
   
   # histogram
   hist(dataset,main=paste("Synthetic Data,",model,"Model"),breaks=50,xlab=xlabel,ylab="Density",freq=F,col="lightgrey")
   
   # theoretical density
   x <- seq(round(min(dataset)),round(max(dataset)),(round(max(dataset))-round(min(dataset)))/500)
   if (model %in% c("LN-LN","rLN-LN")) {
      y <- d.sum.of.mixtures(x,n,p.vector,mu.vector,sigma.vector,logdens=F)   
   }
   else if (model=="EXP-LN") {
      y <- d.sum.of.mixtures(x,n,p.vector,mu.vector,sigma.vector,lambda.vector,logdens=F)   
   }
   
   lines(x,y,col="blue",lwd=3)

   ############
   # estimate #
   ############
   dataset <- as.matrix(dataset,ncol=m)
   result <- stochprof.loop(model=model,dataset=dataset,n=n,TY=TY,fix.mu=F,loops=5,until.convergence=T,print.output=T,show.plots=T,plot.title="Synthetic Data")
   mle <- result$mle
   
   # draw estimated density function
   if (TY>1) {
      p <- mle[1:(TY-1)]
      p <- c(p,1-sum(p))
   }
   else {
      p <- 1
   }
   if (model=="LN-LN") {
      mu <- mle[TY:(length(mle)-1)]
      sigma <- rep(mle[length(mle)],TY)
      z <- d.sum.of.mixtures(x,n,p,mu,sigma,logdens=F)
      
      true.par <- c(p.vector[-length(p.vector)],mu.vector,sigma.vector[1])
   }
   else if (model=="rLN-LN") {
      mu <- mle[TY:((m+1)*TY-1)]
      sigma <- mle[((m+1)*TY):((m+2)*TY-1)]
      z <- d.sum.of.mixtures(x,n,p,mu,sigma,logdens=F)  
      
      true.par <- c(p.vector[-length(p.vector)],mu.vector,sigma.vector)       
   }
   else if (model=="EXP-LN") {
      if (TY>1) {
         mu <- mle[TY:((m+1)*(TY-1))]
         sigma <- rep(mle[(m+1)*(TY-1)+1],TY-1)
      }
      else {
         mu <- NULL
         sigma <- NULL
      }
      lambda <- mle[length(mle)-m+(1:m)]
      z <- d.sum.of.mixtures(x,n,p,mu,sigma,lambda,logdens=F)   
      
      true.par <- c(p.vector[-length(p.vector)],mu.vector,sigma.vector[1],lambda.vector)
   }
   
   hist(dataset,main=paste("Synthetic Data,",model,"Model"),breaks=50,xlab=xlabel,ylab="Density",freq=F,col="lightgrey")   
   lines(x,y,col="blue",lwd=3)   
   lines(x,z,col="darkgreen",lwd=3)
   legend("topright",legend=c("true pdf","estimated pdf"),col=c("blue","darkgreen"),lty=1,lwd=3)
   
   ################
   # print result #
   ################

   names(true.par) <- names(mle)   
   cat("True parameter:\n")
   cat(true.par,"\n\n")

   return(invisible(result))
}
