
calibrate.sigmoid.curve <- function(data , plotit=T , plotspline=F , plotbyds=F) {
  if(class(data)!="list") data <- list(Item1=data)
  NAPPA.output <- lapply(data , NAPPA , output=c("All","Steps"))
  Classes <- lapply(NAPPA.output , function(x) x$Steps$META$Class)
  hks <- lapply(NAPPA.output , function(x) x$HousekeepingFactor)
  geprehks <- lapply(NAPPA.output , function(x) log2(x$Steps$ZEROES[x$Steps$META$Class=="Endogenous",]) )
  coefs <- lapply(1:length(NAPPA.output) , 
                  function(i) t(sapply(1:nrow(geprehks[[i]]) , function(j) coef(summary(lm(geprehks[[i]][j,] ~ hks[[i]])))[2,1:2] )))
  coefs <- lapply(coefs , function(x) { out <- as.data.frame(x); colnames(out) <- c("SLOPE","SE") ; out } )
  for(i in 1:length(coefs)) coefs[[i]]$DS <- i
  for(i in 1:length(coefs)) coefs[[i]]$GeneMeans <- apply(geprehks[[i]] , 1 , mean , na.rm=T)
  for(i in 1:length(coefs)) coefs[[i]]$weights <- 1/coefs[[i]]$SE
    
  coefs.all <- do.call(rbind , coefs)
  
  fit <- nls(SLOPE ~ SSlogis( GeneMeans , Asym = 1, xmid, scal) , 
             start = list(xmid = -3.4, scal = 2) , weights = weights , data=coefs.all)
  parameters <- coef(fit)
  
  if(plotit) {
    plot(coefs.all[,c("GeneMeans","SLOPE")] , pch=16 , col=rgb(0,0,1,0.1) , ylab="Fitted Slopes" , xlab="Gene Averages",axes=F)
    box()
    axis(1)
    axis(2 , at=seq(0,1,0.2) , las=2)
    abline(h=c(0,1))  
    x <- seq(min(coefs.all$GeneMeans) , max(coefs.all$GeneMeans) , length.out=100)
    lines(x , 1/(1+exp((parameters["xmid"]-x)/parameters["scal"])) , lwd=3 , col="RED")
    if(plotspline) {
      lines(x , predict(smooth.spline(x=coefs.all$GeneMeans , y=coefs.all$SLOPE , w=coefs.all$weights , df=4) , x)$y , col="GREEN")      
    }
    if(plotbyds & (length(data)>1)) {
      fitsbyds <- lapply(coefs , function(x) nls(SLOPE ~ SSlogis( GeneMeans , Asym = 1, xmid, scal) , 
                                                           start = list(xmid = -3.4, scal = 2) , weights = weights , data=x ) )
      parametersbyds <- lapply(fitsbyds , coef )
      
      for(i in 1:length(coefs)) {
        lines(x , 1/(1+exp((parametersbyds[[i]]["xmid"]-x)/parametersbyds[[i]]["scal"])) , lwd=1 , lty=2 , col="RED")
      }
      rss.sep <- sum(sapply(fitsbyds , function(x) x$m$deviance()))
      rss.common <- fit$m$deviance()
      Fstat <- ( (rss.common - rss.sep)/(max(coefs.all$DS)-1)) / (rss.sep/(nrow(coefs.all)-max(coefs.all$DS)))
      p <- 1 - pf(Fstat , max(coefs.all$DS)-1 , nrow(coefs.all)-max(coefs.all$DS) )
      legend("bottomright" , paste("p=",round(p,3),sep="") , bty="n")
    }
  }
  
  parameters
}

