lqr = function(y,x,p=0.5,dist = "normal",nu="",gama="",precision = 10^-6,envelope=FALSE,CI=0.95)
{
  par(mfrow=c(1,1))
  if(length(p)==1)
  {
    ## Verify error at parameters specification
    envelope = TRUE
    
    if(dist != "" && dist != "normal" && dist != "t" && dist != "laplace" && dist != "slash" && dist != "cont") stop("The dist values are normal, t, laplace, slash or cont.")
    
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
    if(gama != "" && (gama >= 1 | gama <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
    if(nu != "" && (nu >= 1 | nu <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
    if(nu != "" && (nu >= 100 | nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
    if(nu != "" && nu <= 0 && dist == "slash") stop("nu must be a positive real number.")
    
    if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
    #if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(CI >= 1 | CI <= 0) stop("CI must be a real number in (0,1)")
    if(is.logical(envelope) == FALSE) stop("envelope must be TRUE or FALSE.")
    
    #Running the algorithm
    #out <- suppressWarnings(EM(y,x,p,dist,nu,gama,precision,envelope))
    out <- EM(y,x,p,dist,nu,gama,precision,envelope)
    
    cat('\n')
    cat('--------------------------------------------------------------\n')
    cat('        Quantile Linear Regression using SKD family\n')
    cat('--------------------------------------------------------------\n')
    cat('\n')
    cat("Distribution =",dist,'\n')
    cat("Quantile =",p,'\n')
    #cat("Iterations =",out$iter)
    cat('\n')
    cat('-----------\n')
    cat('Estimates\n')
    cat('-----------\n')
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
    cat('------------------------\n')
    cat('Model selection criteria\n')
    cat('------------------------\n')
    cat('\n')
    critFin <- c(out$loglik, out$AIC, out$BIC, out$HQ)
    critFin <- round(t(as.matrix(critFin)),digits=3)
    dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
    print(critFin)
  
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
  }
  else
  {
    p = sort(p)
    obj.out  = vector("list", length(p))
    
    ## Verify error at parameters specification
    
    if(all(p > 0 && p < 1) == FALSE) stop("p vector must contain real values in (0,1)")
    if(dist != "" && dist != "normal" && dist != "t" && dist != "laplace" && dist != "slash" && dist != "cont") stop("The dist values are normal, t, laplace, slash or cont.")
    ## Verify error at parameters specification
    
    #No data
    if( (length(x) == 0) | (length(y) == 0)) stop("All parameters must be provided.")
    
    #Validating if exists NA's
    if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
    if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
    
    #Validating dims data set
    if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
    if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
    
    #Validating supports
    if(gama != "" && (gama >= 1 | gama <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
    if(nu != "" && (nu >= 1 | nu <= 0) && dist == "cont") stop("nu must be a real number in (0,1)")
    if(nu != "" && (nu >= 100 | nu < 2) && dist == "t") stop("nu must be a positive real number at least 2.")
    if(nu != "" && nu <= 0 && dist == "slash") stop("nu must be a positive real number.")
    
    if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
    #if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(CI >= 1 | CI <= 0) stop("CI must be a real number in (0,1)")
    if(is.logical(envelope) == FALSE) stop("show.convergence must be TRUE or FALSE.")
    
    for(k in 1:length(p))
    {
      #Running the algorithm
      out <- EM(y,x,p[k],dist,nu,gama,precision,envelope)
      
      cat('\n')
      cat('--------------------------------------------------------------\n')
      cat('        Quantile Linear Regression using SKD family\n')
      cat('--------------------------------------------------------------\n')
      cat('\n')
      cat("Distribution =",dist,'\n')
      cat("Quantile =",p[k],'\n')
      #cat("Iterations =",out$iter)
      cat('\n')
      cat('-----------\n')
      cat('Estimates\n')
      cat('-----------\n')
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
      cat('------------------------\n')
      cat('Model selection criteria\n')
      cat('------------------------\n')
      cat('\n')
      critFin <- c(out$loglik, out$AIC, out$BIC, out$HQ)
      critFin <- round(t(as.matrix(critFin)),digits=3)
      dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
      print(critFin)
      
      if(dist == "normal" || dist == "laplace"){
        obj.outk = list(iter = out$iter,criteria = out$criterio,
                       beta = out$theta[1:ncol(as.matrix(x))],
                       sigma= out$theta[ncol(as.matrix(x))+1],
                       SE=out$SE,table = out$table,loglik=out$loglik,
                       AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
        fitted.values = out$fitted.values,residuals=out$residuals)
      }
      
      if(dist == "t" || dist == "slash"){
        obj.outk = list(iter = out$iter,criteria = out$criterio,
                       beta = out$theta[1:ncol(as.matrix(x))],
                       sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
                       SE=out$SE,table = out$table,loglik=out$loglik,
                       AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                       fitted.values = out$fitted.values,residuals=out$residuals)
      }
      
      if(dist == "cont"){
        obj.outk = list(iter = out$iter,criteria = out$criterio,
                       beta = out$theta[1:ncol(as.matrix(x))],
                       sigma= out$theta[ncol(as.matrix(x))+1],nu = out$nu,
                       gamma = out$gamma,
                       SE=out$SE,table = out$table,loglik=out$loglik,
                       AIC=out$AIC,BIC=out$BIC,HQ=out$HQ,
                       fitted.values = out$fitted.values,residuals=out$residuals)
      }
      obj.out[[k]] = obj.outk
    }
    
    par(mfrow=c(1,1))
    d=length(obj.out[[1]]$beta)
    betas = eps = matrix(NA,length(p),d+1)
    
    for (i in 1:length(p))
    {
      j = p[i]
      betas[i,] = c(obj.out[[i]]$beta,obj.out[[i]]$sigma)
      eps[i,] = obj.out[[i]]$SE[1:(d+1)]
    }
    
    LIMSUP = t(betas + qnorm(1-((1-(CI))/2))*eps)
    LIMINF = t(betas - qnorm(1-((1-(CI))/2))*eps)
    labels = list()
    for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
    labels[[d+1]] = bquote(sigma)
    
    par(mar=c(4, 4.5, 1, 0.5))
    op <- par(mfrow=c(ifelse((d+1)%%2==0,(d+1)%/%2,((d+1)%/%2)+1),2),oma=c(0,0,2,0))
    
    for(i in 1:(d+1)){
      
      if(length(p)<4)
      {
        ymin = min(betas[,i],LIMSUP[i,],LIMINF[i,])
        ymax = max(betas[,i],LIMSUP[i,],LIMINF[i,])
        
        plot(p,betas[,i],ylim=c(ymin-2*mean(eps[,i]),ymax+2*mean(eps[,i])),xaxt='n', type='n',xlab='quantiles',ylab=labels[[i]])
        axis(side=1, at=p)
        polygon(c(p,rev(p)),c(LIMSUP[i,],rev(LIMINF[i,])), col = "gray50", border = NA)
        lines(p,betas[,i])
        lines(p,LIMSUP[i,])
        lines(p,LIMINF[i,])
        abline(h=0,lty=2)
      }
      else
      {
        smoothingSpline = smooth.spline(p, betas[,i], spar=0.1)
        smoothingSplineU = smooth.spline(p, betas[,i]+(qnorm(1-((1-(CI))/2)))*eps[,i], spar=0.1)
        smoothingSplineL = smooth.spline(p, betas[,i]-(qnorm(1-((1-(CI))/2)))*eps[,i], spar=0.1)
        plot(p, betas[,i], type='n',xaxt='n',xlab='quantiles',lwd=2,ylim=c(min(smoothingSplineL$y)-2*mean(eps[,i]),max(smoothingSplineU$y)+2*mean(eps[,i])),ylab=labels[[i]])
        axis(side=1, at=p)
        
        #create filled polygon in between the lines
        polygon(c(smoothingSplineL$x,rev(smoothingSplineU$x)),c(smoothingSplineU$y,rev(smoothingSplineL$y)), col = "gray50", border = NA)
        
        #plot lines for high and low range
        lines(p, betas[,i], type='l',lwd=1)
        lines(smoothingSplineU,lwd=1)
        lines(smoothingSplineL,lwd=1)
        abline(h=0,lty=2)
      }
    }
    par(mfrow=c(1,1))
    title("Point estimative and 95% CI for model parameters", outer=TRUE)
    
    class(obj.out)  =  "qr"
    return(obj.out)
  }
}