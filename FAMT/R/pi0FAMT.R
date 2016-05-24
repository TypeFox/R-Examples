pi0FAMT <-
function(model, method = c("smoother", "density"), diagnostic.plot = FALSE){
 convexdensity = function(pvalues) {
  delta <- .00001
  k <- 200
  ny <- 1e-6
  pvalues <- sort(pvalues)
  m <- length(pvalues)
  p.c <- ceiling(100*pvalues)/100
  p.f <- floor(100*pvalues)/100
  t.grid <- (1:100)/100
  x.grid <- (0:100)/100
  t.grid.mat <- matrix(t.grid,ncol=1)
  f.hat <- rep(1,101) #f.hat at the x-grid
  f.hat.p <- rep(1,m) #f.hat at the p-values
  theta.hat <- 0.01*which.max(
  apply(t.grid.mat,1,function(theta) sum((2*(theta-pvalues)*(pvalues<theta)/theta^2))))
  # f.theta.hat at the x-grid
  f.theta.hat <- 2*(theta.hat-x.grid)*(x.grid<theta.hat)/theta.hat^2
  # f.theta.hat at the p-vales
  f.theta.hat.p <- 2*(theta.hat-pvalues)*(pvalues<theta.hat)/theta.hat^2
  i<-1
  j<-0
  z<-1
  thetas <- numeric()
  for(j in 1:k) {
    if (sum((f.hat.p-f.theta.hat.p)/f.hat.p)>0) eps <- 0
    else{
      l <- 0
      u <- 1
      while (abs(u-l)>ny){
        eps <- (l+u)/2
        if (sum(((f.hat.p-f.theta.hat.p)/
        ((1-eps)*f.hat.p+eps*f.theta.hat.p))[f.hat.p>0])<0) l <- eps
        else u <- eps
      }
    }
    #if (theta.hat>0 & eps>0) j<-j+1
    f.hat <- (1-eps)*f.hat + eps*f.theta.hat
    pi.0.hat <- f.hat[101]
    d <- -sum((f.theta.hat.p-f.hat.p)/f.hat.p)
    f.hat.p <- 100*(f.hat[100*p.f+1]-f.hat[100*p.c+1])*(p.c-pvalues)+f.hat[100*p.c+1]
    theta.hat <- 0.01*which.max(apply(t.grid.mat,1,function(theta)
    sum((2*(theta-pvalues)*(pvalues<theta)/theta^2)/f.hat.p)))
    f.theta.hat <- 2*(theta.hat-x.grid)*(x.grid<theta.hat)/theta.hat^2
    f.theta.hat.p <- 2*(theta.hat-pvalues)*(pvalues<theta.hat)/theta.hat^2
    if (sum(f.theta.hat.p/f.hat.p)<sum(1/f.hat.p)){
      theta.hat <- 0
      f.theta.hat <- rep(1,101)
      f.theta.hat.p <- rep(1,m)
    }
    if (sum(thetas==theta.hat)==0){
      thetas[i] <- theta.hat
      thetas <- sort(thetas)
      i <- i + 1
    }
    z <- z+1
  }
  return(f.hat = f.hat)
}
  if (class(model)[1]!="FAMTmodel") stop("Class of model should be FAMTmodel")
   pvalues = model$adjpval
   method <- match.arg(method) # the chosen estimation method
   if (method == "smoother") {
      lambda = seq(0,0.9,0.05)
      # estimates of pi0 for different values of lambda
      vec_pi0 = rep(0, length(lambda))
      for (i in 1:length(lambda)) {
            vec_pi0[i] = sum(pvalues >= lambda[i])/(length(pvalues) * (1 - lambda[i]))
      }
      # spline smoothing
      vec_pi0.spl = smooth.spline(lambda, vec_pi0, df = 3)
      # estimation of pi0 : prediction for the spline model in 
      pi0 = max(0, min(predict(vec_pi0.spl, x = max(lambda))$y,1)) # pi0 belongs to [0,1]
      if (diagnostic.plot) {
            #par(mfrow=c(1,2))
            #get(getOption("device"))()
            plot(lambda, vec_pi0, main = "Smoothing curve employed for pi0 estimation", 
                xlab = "lambda", ylab = "pi0", col = 1, xlim = c(0,1), ylim=c(0,1))
            maxy = max(vec_pi0)
            lines(vec_pi0.spl, col = 3)
            abline(h=pi0, col=2, lty=2)
            axis(2,at=round(pi0,2),cex.axis=0.8, col.axis=2, col=2)
            points(max(lambda), pi0, pch = 4, col = 2, cex=1.2)
      }
   }
   if (method == "density"){
      estim_density = convexdensity(pvalues)
      pi0 = estim_density[length(estim_density)]
   }  
   if (diagnostic.plot) {
        #get(getOption("device"))()
        method_name = paste("method =", method)
        pi0_value = paste("pi0 =", round(pi0, 4))
        dev.new()
        histo = hist(pvalues, freq = FALSE, bre = "FD", main = "Diagnostic Plot: Distribution of p-values and estimated pi0", xlim = c(0, 1), xlab = "p-values")
        y0 = dunif(histo$breaks) * pi0
        lines(histo$breaks, y0, col = 2)
        axis(2,at=round(pi0,2),cex.axis=0.8, col.axis=2, col=2)
        maxy = max(histo$density)
        text(0.5, 5/6 * maxy, method_name)
        text(0.5, 4/6 * maxy, pi0_value, col = 2)        
        if (method == "density"){
            lines(seq(0,1,length=length(estim_density)), estim_density, col=3)
            legend("right",col=3, lty=1, "density estimation")
        }
   }
return(pi0=pi0)
}
