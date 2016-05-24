
####################################
#### Calculate transformed residuals
####################################
cdfU <- function(Y.obs,Y.rep, discrete=FALSE){
  n <- length(Y.obs)
  nc <- ncol(Y.rep)
  Y.tmp <- matrix(rep(Y.obs,times=nc),,nc)
  if(discrete){
    res <- rowMeans(Y.rep <= Y.tmp) - runif(n)*rowMeans(Y.rep == Y.tmp)
  } else res <- rowMeans(Y.rep <= Y.tmp)
  res
}
tranR <- function(Y.obs,Y.rep, discrete=FALSE){
  u <- cdfU(Y.obs,Y.rep, discrete)
  u[u<1e-04] <- 1e-04; u[u>(1-1e-04)] <- 1-1e-04
  qnorm(u)
}
####################################
#### Plot transformed residuals
####################################
plot_etran <- function(e.tran, fig=1:4){
  e.tran <- e.tran[e.tran>-1e5 & e.tran<1e5]
  n <- length(e.tran)
  if(length(fig)==4){
    op <- par(mfrow=c(2,2))
  } else if(length(fig)==3){
    op <- par(mfrow=c(3,1))
  } else op <- par(mfrow=c(1,length(fig)))
  
  if(any(fig==1)) {
    plot(0,0, xlim=c(0,n), ylim=c(-3,3), type="n", xlab="observation", ylab="etran", main="Transformed Residuals")
    points(e.tran)
  }
  if(any(fig==2)){
    qqnorm(e.tran); qqline(e.tran,col=2)
  }
  if(any(fig==3)){
    e.den <- density(e.tran)
    plot(x <- seq(-5,5,by=0.1),dnorm(x),col=2, type="n", xlab="x", ylab="f(x)", main="Densities")
    lines(e.den$x, e.den$y)
    lines(x, dnorm(x),col=2, type="l", lty=2, lwd=3)
  }
  if(any(fig==4)){
    if (requireNamespace("reldist", quietly = TRUE)) {
      rd <- reldist::reldist( y=e.tran, yo=qnorm((1:n+0.5)/(n+1)), ci=FALSE, graph=FALSE)
      plot(0,0, xlim=c(0,1), ylim=c(0,2), type="n", xlab="u", ylab="r(u)", main="Relative Densities")
      lines(rd$x, rd$y)
      abline(h=1, col=2, lty=2, lwd=3)
    } else{
      par(op)
      stop("Please install and load {reldist} first before using this function!")
    }
  }
  par(op)
}
####################################
#### Calculate distances from tranR
####################################
e2dist <- function(e.tran)
{
    x <- e.tran
    x <- x[x>-1e5 & x<1e5]
    res <- NULL
    if (requireNamespace("distr", quietly = TRUE) && requireNamespace("distrEx", quietly = TRUE)) {
      res <- c(distrEx::HellingerDist(distr::Norm(), x), 
               distrEx::HellingerDist(x, distr::Norm(), asis.smooth.discretize = "smooth"), 
               distrEx::KolmogorovDist(distr::Norm(),x) )
      names(res) <- c("Hellinger.discre","Hellinger.smooth","Kolmogorov")
    } else{
      stop("Please install and load {distr} and {distrEx} first before using this function!")
    }    
    res
  }
####################################
#### Generate baseline samples
####################################
baseline.dist <- function(n, iter){
res <- matrix(0, iter, 3)
if (requireNamespace("distr", quietly = TRUE) && requireNamespace("distrEx", quietly = TRUE)) {
  for(i in 1:iter){
    x <- rnorm(n)
    res[i,] <- c(distrEx::HellingerDist(x, distr::Norm()), 
                 distrEx::HellingerDist(x, distr::Norm(), asis.smooth.discretize = "smooth"),
                 distrEx::KolmogorovDist(x, distr::Norm())
    )
  }
} else{
  stop("Please install and load {distr} and {distrEx} first before using this function!")
}
colnames(res) <- c("Hellinger.discre","Hellinger.smooth","Kolmogorov")
res
}

# baseline.parallel <- function(n, iter, n.cores = getOption("cores")){
# res0 <- mclapply(1:iter, function(t){
#         x <- rnorm(n)
#         c(HellingerDist(x, Norm()),HellingerDist(x, Norm(), asis.smooth.discretize = "smooth"),KolmogorovDist(x, Norm()) )
#       }, mc.cores = n.cores
#     )
# res <- matrix(unlist(res0),,3,byrow=TRUE)
# colnames(res) <- c("Hellinger.discre","Hellinger.smooth","Kolmogorov")
# res
# }
####################################
#### Plot baseline of distances
####################################
plot_baseline <- function(d.samples, dist.name){
val.crit <- quantile(d.samples, prob=c(0.05,0.5,0.95))
plot(density(d.samples),main=paste("Baseline of ",dist.name," Distance",sep=""))
abline(v=val.crit,lty=2,col=c(3,2,3))
#text(val.crit,0,paste(c("2.5%","50%","97.5%"),round(val.crit,4),sep="="),col=c(3,2,3),cex=0.8)
legend("topright",paste(c("  5%","50%","95%"),round(val.crit,3),sep="="),lty=2,col=c(3,2,3))
}
####################################
#### Calculate one-side p-value
####################################
pOne <- function(d.obs, d.base){
  if(length(d.obs)==1){
  p <- mean(d.base<=d.obs)
  p <- ifelse(p<=0.5, p, 1-p)
  } else{
      tmp <- rbind(d.obs, d.base)
      p <- apply(tmp, 2, function(t){
              pp <- mean(t[-1]<=t[1]) 
              ifelse(pp<=0.5, pp, 1-pp)
              } 
            )
    }
  p
  }
####################################
#### END
####################################