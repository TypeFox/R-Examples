grid_example = function()
{
  xq = seq(5,95,5)
  
  d=3
  betas = eps = matrix(NA,length(xq),d+1)
  
  for (i in 1:length(xq))
  {
    j = xq[i]
    betas[i,] = c(t(get(paste("pitui",j,sep=""))$res$beta),get(paste("pitui",j,sep=""))$res$sigmae)
    eps[i,] = get(paste("pitui",j,sep=""))$res$EP[1:(d+1)]
  }
  
  labels = list()
  for(i in 1:d){labels[[i]] = bquote(beta[.(i)])}
  labels[[d+1]] = bquote(sigma)
  
  par(mar=c(4, 4.5, 1, 0.5))
  op <- par(mfrow=c(ifelse((d+1)%%2==0,(d+1)%/%2,((d+1)%/%2)+1),2),oma=c(0,0,2,0))
  
  
  for(i in 1:4){
    smoothingSpline = smooth.spline(xq/100, betas[,i], spar=0.1)
    #plot(xq/100, betas[,i])
    #lines(smoothingSpline)
    
    smoothingSplineU = smooth.spline(xq/100, betas[,i]+(1.96)*eps[,i], spar=0.1)
    #plot(xq/100, betas[,i]+(1.96)*eps[,i])
    #lines(smoothingSplineU)
    
    smoothingSplineL = smooth.spline(xq/100, betas[,i]-(1.96)*eps[,i], spar=0.1)
    #plot(xq/100, betas[,i]-(1.96)*eps[,i])
    #lines(smoothingSplineL)
    
    
    #Set up plot but don't plot any points (type = 'n')
    plot(xq/100, betas[,i], type='n',xaxt='n',xlab='quantiles',lwd=2,ylim=c(min(smoothingSplineL$y)-2*mean(eps[,i]),max(smoothingSplineU$y)+2*mean(eps[,i])),ylab=labels[[i]])
    axis(side=1, at=xq[2:20]/100)
    #create filled polygon in between the lines
    polygon(c(smoothingSplineL$x,rev(smoothingSplineU$x)),c(smoothingSplineU$y,rev(smoothingSplineL$y)), col = "grey", border = NA)
    
    #plot lines for high and low range
    lines(xq/100, betas[,i], type='l',lwd=1)
    lines(smoothingSplineU,lwd=1)
    lines(smoothingSplineL,lwd=1)
    abline(h=0,lty=2)
  }
  title("Point estimative and 95% CI for model parameters", outer=TRUE)
  
}