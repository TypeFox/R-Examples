best.lqr = function(y,x,p=0.5,precision = 10^-6,criterion = "AIC")
{
  if(length(p)>1) stop("The function best.lqr is only available for one quantile.")
  ## Verify error at parameters specification
  envelope = TRUE
  
  #No data
  if( (length(x) == 0) | (length(y) == 0)) stop("All parameters must be provided.")
  
  #Validating if exists NA's
  if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
  if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
  
  #Validating dims data set
  if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
  if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
  
  #Validating supports
  if(p >= 1 | p <= 0) stop("p must be a real number in (0,1)")
  if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
  if(criterion != "" && criterion != "AIC" && criterion != "BIC" && criterion != "HQ" && criterion != "loglik") stop("The criterion must be AIC, BIC, HQ or loglik.")
  
  #Running the algorithm
  #out <- suppressWarnings(EM(y,x,p,dist,nu,gama,precision,envelope))
  
  vdist = c("normal","t","laplace","slash","cont")
  vDIST = c("Normal","Student-t","Laplace","Slash","Cont. Normal")
  obj.out  = vector("list", 5)
  
  #TEST
  p = 0.5
  precision = 10^-6
  #######################
  #oma margenes totales
  
  par(mfrow=c(2,3),oma=c(0,1,1,0),mar=c(5, 4, 4, 2.5))
  for(k in 1:5)
  {
    obj.out[[k]] <- EM2(y = y,x = x,p = p,dist = vdist[k],precision = precision,envelope=TRUE)
  }
  plot.new()
  #mtext("Histogram of residuals and fitted densities", side = 3, line = 1, outer = TRUE,cex=1.3)
  
  
  #HISTOGRAM AND FITTED DENSITIES
  
  resall = unlist(lapply(X = obj.out,FUN = residuals))
  par(mfrow=c(2,3),oma=c(0,1,0,1),mar=c(5, 4, 5, 0.5))
  seqq = seq(from=min(resall),to = max(resall),length.out = 1000)
  up = dSKD(y = 0,mu = 0,sigma = obj.out[[3]]$theta[dim(x)[2]+1],p = p,dist = vdist[3])
  
  folginha = 0.1
  for(k in 1:5)
  {
    hist(x = obj.out[[k]]$residuals,freq=FALSE,breaks=sqrt(length(y)),xlab = "residuals",
         main = vDIST[k],cex.main=1.5,ylim=c(0,(1+folginha)*up))
    dens = dSKD(y = seqq,mu = 0,sigma = obj.out[[k]]$theta[dim(x)[2]+1],p = p,dist = vdist[k],
                nu=obj.out[[k]]$nu,gama = obj.out[[k]]$gamma)
    lines(seqq,dens,lwd=1.5,col="blue")
  }
  plot.new()
  #mtext("Histogram of residuals and fitted densities", side = 3, line = 1, outer = TRUE,cex=1.3)
  
  RES = matrix(data = NA,nrow = 4,ncol = 5)
  for(k in 1:5)
  {
    RES[1,k] = (obj.out[[k]])$AIC
    RES[2,k] = (obj.out[[k]])$BIC
    RES[3,k] = (obj.out[[k]])$HQ
    RES[4,k] = (obj.out[[k]])$loglik
  }
  colnames(RES) = c(vDIST[1:4],"C. Normal")
  rownames(RES) = c("AIC","BIC","HQ","loglik")
  
  if(criterion == "AIC")
  {
    index = which(RES[1,] == min(RES[1,]))
  }
  if(criterion == "BIC")
  {
    index = which(RES[2,] == min(RES[2,])) 
  }
  if(criterion == "HQ")
  {
    index = which(RES[3,] == min(RES[3,])) 
  }
  if(criterion == "loglik")
  {
    index = which(RES[4,] == max(RES[4,])) 
  }
  
  out  = obj.out[[index]]
  dist = vdist[index] 
  cat('\n')
  cat('--------------------------------------------------------------\n')
  cat('        Quantile Linear Regression using SKD family\n')
  cat('--------------------------------------------------------------\n')
  cat('\n')
  cat("Criterion =",criterion,'\n')
  cat("Best fit =",vDIST[index],'\n')
  cat("Quantile =",p,'\n')
  cat('\n')
  cat('--------------------------------\n')
  cat('Model Likelihood-Based criterion\n')
  cat('--------------------------------\n')
  cat('\n')
  print(RES)
  #cat("Iterations =",out$iter)
  cat('\n')
  cat('---------\n')
  cat('Estimates\n')
  cat('---------\n')
  cat('\n')
  print(out$table)
  cat('---\n')
  cat('Signif. codes:  0 "***" 0.001 "**" 0.01 "*" 0.05 "." 0.1 " " 1\n')
  cat('\n')
  cat('sigma =',round(out$theta[ncol(as.matrix(x))+1],5),'\n')
  if(dist == "normal" || dist == "laplace"){
    cat('\n')
  }
  if(dist == "t" || dist == "slash"){
    cat('nu    =',out$nu,'\n')
    cat('\n')
  }
  if(dist == "cont"){
    cat('nu    =',out$nu,'\n')
    cat('gamma =',out$gamma,'\n')
    cat('\n')
  }
  
  if(dist == "normal" || dist == "laplace"){
    obj.out = list(iter = out$iter,criteria = out$criterio,
                   beta = out$theta[1:ncol(as.matrix(x))],
                   sigma= out$theta[ncol(as.matrix(x))+1],
                   SE=out$SE,table = out$table,loglik=out$loglik,
                   AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                   fitted.values = out$fitted.values,residuals=out$residuals)
  }
  
  if(dist == "t" || dist == "slash"){
    obj.out = list(iter = out$iter,criteria = out$criterio,
                   beta = out$theta[1:ncol(as.matrix(x))],
                   sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
                   SE=out$SE,table = out$table,loglik=out$loglik,
                   AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                   fitted.values = out$fitted.values,residuals=out$residuals)
  }
  
  if(dist == "cont"){
    obj.out = list(iter = out$iter,criteria = out$criterio,
                   beta = out$theta[1:ncol(as.matrix(x))],
                   sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
                   gamma = out$gamma,
                   SE=out$SE,table = out$table,loglik=out$loglik,
                   AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                   fitted.values = out$fitted.values,residuals=out$residuals)
  }
  class(obj.out)  =  "qr"
  return(obj.out)  
  
  # else
  # {
  #   p = sort(p)
  #   obj.out  = vector("list", length(p))
  #   
  #   ## Verify error at parameters specification
  #   
  #   if(all(p > 0 && p < 1) == FALSE) stop("p vector must contain real values in (0,1)")
  #   
  #   ## Verify error at parameters specification
  #   
  #   #No data
  #   if( (length(x) == 0) | (length(y) == 0)) stop("All parameters must be provided.")
  #   
  #   #Validating if exists NA's
  #   if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
  #   if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
  #   
  #   #Validating dims data set
  #   if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
  #   if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
  #   
  #   #Validating supports
  #   if(gama != "" && (gama >= 1 | gama <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  #   if(nu != "" && (nu >= 1 | nu <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
  #   if(nu != "" && (nu >= 100 | nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
  #   if(nu != "" && nu <= 0 && dist == "slash") stop("nu must be a positive real number.")
  #   
  #   if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
  #   #if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
  #   if(CI >= 1 | CI <= 0) stop("CI must be a real number in (0,1)")
  #   if(is.logical(envelope) == FALSE) stop("show.convergence must be TRUE or FALSE.")
  #   
  #   for(k in 1:length(p))
  #   {
  #     #Running the algorithm
  #     out <- EM(y,x,p[k],dist,nu,gama,precision,envelope)
  #     
  #     cat('\n')
  #     cat('--------------------------------------------------------------\n')
  #     cat('        Quantile Linear Regression using SKD family\n')
  #     cat('--------------------------------------------------------------\n')
  #     cat('\n')
  #     cat("Distribution =",dist,'\n')
  #     cat("Quantile =",p[k],'\n')
  #     #cat("Iterations =",out$iter)
  #     cat('\n')
  #     cat('-----------\n')
  #     cat('Estimates\n')
  #     cat('-----------\n')
  #     cat('\n')
  #     print(out$table)
  #     cat('---\n')
  #     cat('Signif. codes:  0 "***" 0.001 "**" 0.01 "*" 0.05 "." 0.1 " " 1\n')
  #     cat('\n')
  #     cat('sigma =',round(out$theta[ncol(as.matrix(x))+1],5),'\n')
  #     if(dist == "normal" || dist == "laplace"){
  #       cat('\n')
  #     }
  #     if(dist == "t" || dist == "slash"){
  #       cat('nu    =',out$nu,'\n')
  #       cat('\n')
  #     }
  #     if(dist == "cont"){
  #       cat('nu    =',out$nu,'\n')
  #       cat('gamma =',out$gamma,'\n')
  #       cat('\n')
  #     }
  #     cat('------------------------\n')
  #     cat('Model selection criteria\n')
  #     cat('------------------------\n')
  #     cat('\n')
  #     critFin <- c(out$loglik, out$AIC, out$BIC, out$HQ)
  #     critFin <- round(t(as.matrix(critFin)),digits=3)
  #     dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
  #     print(critFin)
  #     
  #     if(dist == "normal" || dist == "laplace"){
  #       obj.outk = list(iter = out$iter,criteria = out$criterio,
  #                       beta = out$theta[1:ncol(as.matrix(x))],
  #                       sigma= out$theta[ncol(as.matrix(x))+1],
  #                       SE=out$SE,table = out$table,loglik=out$loglik,
  #                       AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
  #                       fitted.values = out$fitted.values,residuals=out$residuals)
  #     }
  #     
  #     if(dist == "t" || dist == "slash"){
  #       obj.outk = list(iter = out$iter,criteria = out$criterio,
  #                       beta = out$theta[1:ncol(as.matrix(x))],
  #                       sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
  #                       SE=out$SE,table = out$table,loglik=out$loglik,
  #                       AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
  #                       fitted.values = out$fitted.values,residuals=out$residuals)
  #     }
  #     
  #     if(dist == "cont"){
  #       obj.outk = list(iter = out$iter,criteria = out$criterio,
  #                       beta = out$theta[1:ncol(as.matrix(x))],
  #                       sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
  #                       gamma = out$gamma,
  #                       SE=out$SE,table = out$table,loglik=out$loglik,
  #                       AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
  #                       fitted.values = out$fitted.values,residuals=out$residuals)
  #     }
  #     obj.out[[k]] = obj.outk
  #   }
  #   
  #   par(mfrow=c(1,1))
  #   d=length(obj.out[[1]]$beta)
  #   betas = eps = matrix(NA,length(p),d+1)
  #   
  #   for (i in 1:length(p))
  #   {
  #     j = p[i]
  #     betas[i,] = c(obj.out[[i]]$beta,obj.out[[i]]$sigma)
  #     eps[i,] = obj.out[[i]]$SE[1:(d+1)]
  #   }
  #   
  #   LIMSUP = t(betas + qnorm(1-((1-(CI))/2))*eps)
  #   LIMINF = t(betas - qnorm(1-((1-(CI))/2))*eps)
  #   labels = list()
  #   for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
  #   labels[[d+1]] = bquote(sigma)
  #   
  #   par(mar=c(4, 4.5, 1, 0.5))
  #   op <- par(mfrow=c(ifelse((d+1)%%2==0,(d+1)%/%2,((d+1)%/%2)+1),2),oma=c(0,0,2,0))
  #   
  #   for(i in 1:(d+1)){
  #     
  #     if(length(p)<4)
  #     {
  #       ymin = min(betas[,i],LIMSUP[i,],LIMINF[i,])
  #       ymax = max(betas[,i],LIMSUP[i,],LIMINF[i,])
  #       
  #       plot(p,betas[,i],ylim=c(ymin-2*mean(eps[,i]),ymax+2*mean(eps[,i])),xaxt='n', type='n',xlab='quantiles',ylab=labels[[i]])
  #       axis(side=1, at=p)
  #       polygon(c(p,rev(p)),c(LIMSUP[i,],rev(LIMINF[i,])), col = "gray50", border = NA)
  #       lines(p,betas[,i])
  #       lines(p,LIMSUP[i,])
  #       lines(p,LIMINF[i,])
  #       abline(h=0,lty=2)
  #     }
  #     else
  #     {
  #       smoothingSpline = smooth.spline(p, betas[,i], spar=0.1)
  #       smoothingSplineU = smooth.spline(p, betas[,i]+(qnorm(1-((1-(CI))/2)))*eps[,i], spar=0.1)
  #       smoothingSplineL = smooth.spline(p, betas[,i]-(qnorm(1-((1-(CI))/2)))*eps[,i], spar=0.1)
  #       plot(p, betas[,i], type='n',xaxt='n',xlab='quantiles',lwd=2,ylim=c(min(smoothingSplineL$y)-2*mean(eps[,i]),max(smoothingSplineU$y)+2*mean(eps[,i])),ylab=labels[[i]])
  #       axis(side=1, at=p)
  #       
  #       #create filled polygon in between the lines
  #       polygon(c(smoothingSplineL$x,rev(smoothingSplineU$x)),c(smoothingSplineU$y,rev(smoothingSplineL$y)), col = "gray50", border = NA)
  #       
  #       #plot lines for high and low range
  #       lines(p, betas[,i], type='l',lwd=1)
  #       lines(smoothingSplineU,lwd=1)
  #       lines(smoothingSplineL,lwd=1)
  #       abline(h=0,lty=2)
  #     }
  #   }
  #   par(mfrow=c(1,1))
  #   title("Point estimative and 95% CI for model parameters", outer=TRUE)
  #   
  #   class(obj.out)  =  "qr"
  #   return(obj.out)
  # }
}