QRLMM = function(y,x,z,nj,p=0.5,precision=0.0001,MaxIter=300,M=10,cp=0.25,beta=NA,sigma=NA,Psi=NA,show.convergence=TRUE,CI=95)
{
  if(length(p)==1)
  {
    ## Verify error at parameters specification
    
    #No data
    if( (length(x) == 0) | (length(y) == 0) | (length(z) == 0)) stop("All parameters must be provided.")
    
    #Validating if exists NA's
    if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
    if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
    if(sum(y[is.na(z)==TRUE]) > 0) stop("There are some NA's values in z")
    
    #Validating dims data set
    if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
    if( length(y) != sum(nj) ) stop("nj does not match with  the provided data. (length(y) != sum(nj))")
    if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
    if( length(y) != nrow(as.matrix(z)) ) stop("z variable does not have the same number of lines than y")
    
    #Validating supports
    if(p >= 1 | p <= 0) stop("p must be a real number in (0,1)")
    if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
    if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(M <= 0 |M%%1!=0) stop("M must be a positive integer value >= 10")
    if(cp > 1 | cp < 0) stop("cp must be a real number in [0,1]")
    if(is.logical(show.convergence) == FALSE) stop("show.convergence must be TRUE or FALSE.")
    
    #Matrix column labels
    namesz <- ('z1')
    if(ncol(as.matrix(z))>1){
      for(i in 2:ncol(as.matrix(z))){namesz <- c(namesz, paste("z",i,sep=""))}
    }
    
    #FALTA VALIDAR VALORES INICIAIS
    
    #No data
    if(is.na(beta) == FALSE)
    {if(length(beta) != ncol(as.matrix(x))) stop("beta must be a vector of dimension equal to columns of x")}
    
    if(is.na(sigma) == FALSE)
    {if(sigma <= 0) stop("sigma must be a positive number")}
    
    #Load required libraries
    
    if(is.na(Psi) == FALSE)
    {  
      if( (ncol(as.matrix(D)) != ncol(as.matrix(z))) | ((nrow(as.matrix(D)) != ncol(as.matrix(z)))) ) stop("D must be a square matrix of dims equal to the columns of z")
      if(is.positive.definite(D)==FALSE) stop("D must be a square symmetrical real posite definite matrix.") 
    }
    
    #intial values
    if(is.na(beta) == TRUE){beta   = suppressWarnings(rq(y ~ -1 + x,tau = p))$coefficients}
    if(is.na(sigma) == TRUE){dif   = y - x%*%beta;
                             sigma = (1/length(y))*(sum(p*dif[dif>0]) - sum((1-p)*dif[dif<0]))}
    if(is.na(Psi) == TRUE){Psi = diag(ncol(as.matrix(z)))}
    
    #Running the algorithm
    out <- suppressWarnings(QSAEM_COM_7(y,x,z,nj,p,precision,MaxIter,M,pc=cp,beta=beta,sigmae=sigma,D=Psi))
    
    cat('\n')
    cat('---------------------------------------------------\n')
    cat('Quantile Regression for Linear Mixed Model\n')
    cat('---------------------------------------------------\n')
    cat('\n')
    cat("Quantile =",p)
    cat('\n')
    cat("Subjects =",length(nj),";",'Observations =',sum(nj),
        ifelse(sum(nj==nj[1])==length(nj),'; Balanced =',""),
        ifelse(sum(nj==nj[1])==length(nj),nj[1],""))
    cat('\n')
    cat('\n')
    cat('-----------\n')
    cat('Estimates\n')
    cat('-----------\n')
    cat('\n')
    cat('- Fixed effects \n')
    cat('\n')
    print(round(out$res$table,5))
    cat('\n')
    cat('sigma =',round(out$res$sigmae,5),'\n')
    cat('\n')
    cat('Random effects Variance-Covariance Matrix matrix \n')
    dimnames(out$res$D) <- list(namesz,namesz)
    print(round(out$res$D,5))
    cat('\n')
    cat('------------------------\n')
    cat('Model selection criteria\n')
    cat('------------------------\n')
    cat('\n')
    critFin <- c(out$res$loglik, out$res$AIC, out$res$BIC, out$res$HQ)
    critFin <- round(t(as.matrix(critFin)),digits=3)
    dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
    print(critFin)
    cat('\n')
    cat('-------\n')
    cat('Details\n')
    cat('-------\n')
    cat('\n')
    cat("Convergence reached? =",(out$res$iter < MaxIter))
    cat('\n')
    cat('Iterations =',out$res$iter,"/",MaxIter)
    cat('\n')
    cat('Criteria =',round(out$res$criterio,5))
    cat('\n')
    cat('MC sample =',M)
    cat('\n')  
    cat('Cut point =',cp)
    cat('\n')
    cat("Processing time =",out$res$time,units(out$res$time))
    
    if(show.convergence=="TRUE")
    {
      cpl = cp*MaxIter
      d = dim(x)[2]
      q = dim(z)[2]
      ndiag  = (q*(1+q)/2)
      npar   = d+1+ndiag
      
      labels = list()
      for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
      labels[[d+1]] = bquote(sigma)
      for(i in 1:ndiag){labels[[i+d+1]] = bquote(psi[.(i)])}
      
      par(mar=c(4, 4.5, 1, 0.5))
      op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3))
      
      for(i in 1:npar)
      {
        plot.ts(out$conv$teta[i,],xlab="Iteration",ylab=labels[[i]])
        abline(v=cpl,lty=2)
      }
    }
    
    res     = list(iter = out$res$iter,criteria = out$res$criterio,beta = out$res$beta,sigma= out$res$sigmae,Psi = out$res$D,SE=out$res$EP,table = out$res$table,loglik=out$res$loglik,AIC=out$res$AIC,BIC=out$res$BIC,HQ=out$res$HQ,time = out$res$time)
    obj.out = list(conv=out$conv,res = res)
    class(obj.out)  =  "QRLMM"
    return(obj.out)  
  }
  else
  {
    p = sort(p)
    obj.out  = vector("list", length(p))
    
    ## Verify error at parameters specification
    
    #No data
    if( (length(x) == 0) | (length(y) == 0) | (length(z) == 0)) stop("All parameters must be provided.")
    
    #Validating if exists NA's
    if(sum(y[is.na(y)==TRUE]) > 0) stop("There are some NA's values in y")
    if(sum(y[is.na(x)==TRUE]) > 0) stop("There are some NA's values in x")
    if(sum(y[is.na(z)==TRUE]) > 0) stop("There are some NA's values in z")
    
    #Validating dims data set
    if(ncol(as.matrix(y)) > 1) stop("y must have just one column")
    if( length(y) != sum(nj) ) stop("nj does not match with  the provided data. (length(y) != sum(nj))")
    if( length(y) != nrow(as.matrix(x)) ) stop("x variable does not have the same number of lines than y")
    if( length(y) != nrow(as.matrix(z)) ) stop("z variable does not have the same number of lines than y")
    
    #Validating supports
    if(all(p > 0 && p < 1) == FALSE) stop("p vector must contain real values in (0,1)")
    if(precision <= 0) stop("precision must be a positive value (suggested to be small)")
    if(MaxIter <= 0 |MaxIter%%1!=0) stop("MaxIter must be a positive integer value")
    if(M <= 0 |M%%1!=0) stop("M must be a positive integer value >= 10")
    if(cp >= 1 | cp <= 0) stop("cp must be a real number in [0,1]")
    if(is.logical(show.convergence) == FALSE) stop("show.convergence must be TRUE or FALSE.")
    
    #Matrix column labels
    namesz <- ('z1')
    if(ncol(as.matrix(z))>1){
      for(i in 2:ncol(as.matrix(z))){namesz <- c(namesz, paste("z",i,sep=""))}
    }
    
    #Load required libraries
    pb2 = tkProgressBar(title = "QRLMM for several quantiles", min = 0,max = length(p), width = 300)
    
    for(k in 1:length(p))
    {
      setTkProgressBar(pb2, k-1, label=paste("Running quantile ",p[k],"   -   ",k-1,"/",length(p),"   -   ",round((k-1)/length(p)*100,0),"% done",sep = ""))
      
      beta   = suppressWarnings(rq(y ~ -1 + x,tau = p[k]))$coefficients
      dif   = y - x%*%beta
      sigma = (1/length(y))*(sum(p[k]*dif[dif>0]) - sum((1-p[k])*dif[dif<0]))
      Psi = diag(ncol(as.matrix(z)))
      
      #Running the algorithm
      out <- suppressWarnings(QSAEM_COM_7(y,x,z,nj,p[k],precision,MaxIter,M,pc=cp,beta=beta,sigmae=sigma,D=Psi))
      
      cat('\n')
      cat('---------------------------------------------------\n')
      cat('Quantile Regression for Linear Mixed Model\n')
      cat('---------------------------------------------------\n')
      cat('\n')
      cat("Quantile =",p)
      cat('\n')
      cat("Subjects =",length(nj),";",'Observations =',sum(nj),
          ifelse(sum(nj==nj[1])==length(nj),'; Balanced =',""),
          ifelse(sum(nj==nj[1])==length(nj),nj[1],""))
      cat('\n')
      cat('\n')
      cat('-----------\n')
      cat('Estimates\n')
      cat('-----------\n')
      cat('\n')
      cat('- Fixed effects \n')
      cat('\n')
      print(round(out$res$table,5))
      cat('\n')
      cat('sigma =',round(out$res$sigmae,5),'\n')
      cat('\n')
      cat('Random effects Variance-Covariance Matrix matrix \n')
      dimnames(out$res$D) <- list(namesz,namesz)
      print(round(out$res$D,5))
      cat('\n')
      cat('------------------------\n')
      cat('Model selection criteria\n')
      cat('------------------------\n')
      cat('\n')
      critFin <- c(out$res$loglik, out$res$AIC, out$res$BIC, out$res$HQ)
      critFin <- round(t(as.matrix(critFin)),digits=3)
      dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQ"))
      print(critFin)
      cat('\n')
      cat('-------\n')
      cat('Details\n')
      cat('-------\n')
      cat('\n')
      cat("Convergence reached? =",(out$res$iter < MaxIter))
      cat('\n')
      cat('Iterations =',out$res$iter,"/",MaxIter)
      cat('\n')
      cat('Criteria =',round(out$res$criterio,5))
      cat('\n')
      cat('MC sample =',M)
      cat('\n')  
      cat('Cut point =',cp)
      cat('\n')
      cat("Processing time =",out$res$time,units(out$res$time))
      
      res      = list(iter = out$res$iter,criteria = out$res$criterio,beta = out$res$beta,sigma= out$res$sigmae,Psi = out$res$D,SE=out$res$EP,table = out$res$table,loglik=out$res$loglik,AIC=out$res$AIC,BIC=out$res$BIC,HQ=out$res$HQ,time = out$res$time)
      obj.outk = list(conv=out$conv,res = res)
      obj.out[[k]] = obj.outk
    }
    close(pb2)
    
    par(mfrow=c(1,1))
    d=length(obj.out[[1]]$res$beta)
    betas = eps = matrix(NA,length(p),d+1)
    
    for (i in 1:length(p))
    {
      j = p[i]
      betas[i,] = rbind(obj.out[[i]]$res$beta,obj.out[[i]]$res$sigma)
      eps[i,] = obj.out[[i]]$res$SE[1:(d+1)]
    }
    
    LIMSUP = t(betas + qnorm(1-((1-(CI/100))/2))*eps)
    LIMINF = t(betas - qnorm(1-((1-(CI/100))/2))*eps)
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
      }
      else
      {
        smoothingSpline = smooth.spline(p, betas[,i], spar=0.1)
        smoothingSplineU = smooth.spline(p, betas[,i]+(qnorm(1-((1-(CI/100))/2)))*eps[,i], spar=0.1)
        smoothingSplineL = smooth.spline(p, betas[,i]-(qnorm(1-((1-(CI/100))/2)))*eps[,i], spar=0.1)
        plot(p, betas[,i], type='n',xaxt='n',xlab='quantiles',lwd=2,ylim=c(min(smoothingSplineL$y)-2*mean(eps[,i]),max(smoothingSplineU$y)+2*mean(eps[,i])),ylab=labels[[i]])
        axis(side=1, at=p)
        
        #create filled polygon in between the lines
        polygon(c(smoothingSplineL$x,rev(smoothingSplineU$x)),c(smoothingSplineU$y,rev(smoothingSplineL$y)), col = "gray50", border = NA)
        
        #plot lines for high and low range
        lines(p, betas[,i], type='l',lwd=1)
        lines(smoothingSplineU,lwd=1)
        lines(smoothingSplineL,lwd=1)
      }
    }
    title("Point estimative and 95% CI for model parameters", outer=TRUE)
    
    if(show.convergence=="TRUE")
    {
      cpl = cp*MaxIter
      d = dim(x)[2]
      q = dim(z)[2]
      ndiag  = (q*(1+q)/2)
      npar   = d+1+ndiag
      
      labels = list()
      for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
      labels[[d+1]] = bquote(sigma)
      for(i in 1:ndiag){labels[[i+d+1]] = bquote(psi[.(i)])}
      
      for(k in 1:length(p))
      {
        cat('\n')
        cat ("Press [ENTER] to see next plot:")
        line <- readline()
        
        par(mar=c(4, 4.5, 1, 0.5))
        op <- par(mfrow=c(ifelse(npar%%3==0,npar%/%3,(npar%/%3)+1),3),oma=c(0,0,2,0))
        
        for(i in 1:npar)
        {
          plot.ts(obj.out[[k]]$conv$teta[i,],xlab="Iteration",ylab=labels[[i]])
          abline(v=cpl,lty=2)
        }
        title(paste("Convergence plots for quantile",p[k]), outer=TRUE)
      }
    }
    class(obj.out)  =  "QRLMM"
    return(obj.out)
  }
}