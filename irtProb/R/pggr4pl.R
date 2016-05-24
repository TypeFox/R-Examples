`pggr4pl` <-
function(x=ggr4pl(rep=1), rep=1, n=dim(x)[2], N=dim(x)[1], theta=rep(0,N), a=rep(1,n), b=rep(0,n), c=rep(0,n), d=rep(1,n), log.p=FALSE, TCC=FALSE) {
 N <- N # Unused argument that would have to be retired
 p           <- NULL
 for (i in 1:n) p <- cbind(p, p4pl(theta=theta,a=a[i],b=b[i],c=c[i],d=d[i]))
 q           <- 1-p
 prob        <- p^x * q^(1-x)
  # La fonction rep() permet de s'assurer du bon ordre de theta
 prob        <- data.frame(theta=rep(theta, each=rep), total=rowSums(x), prob=apply(prob,1,prod))
 if (TCC==TRUE) {
  figure   <- NULL
  simTheta <- seq(-4,4,length=50)
  for (i in 1:N) {
   #simTheta <- seq(-4+theta[i],4+theta[i],length=100)
   simProb    <- temp <- NULL
   true.theta <- prob$theta[i]
   for (j in simTheta) {
    temp    <- cbind(true.theta, pggr4pl(x[i,],n=n,N=1, theta=j, a=a,b=b,c=c,d=d))
    simProb <- rbind(simProb, temp)
    }
   simProb <- cbind(N=rep(i,length(simTheta)), simProb)
   figure  <- rbind(figure, simProb)
   }
  # Voir Sarkar (2008, p. 71-74) pour l'utilisation de subscripts pour indicer les valeurs de la conditioning variable
   tcc   <- xyplot(prob ~  theta | factor(N), data=figure, stats=figure$true.theta,
                  xlab="Theta", ylab="Probability", par.strip.text=list(cex=0.80),
                  panel = function(x,y,stats, ..., subscripts) {
                   panel.xyplot(x,y,                       lty="solid",  type="l")
                   panel.abline(v=sapply(stats[subscripts], mean), lty="dashed", col="red")  ##???
                   }
                  )
  }
  if (log.p==TRUE) prob$prob <- log(prob$prob)
  if (TCC ==TRUE)  return(list(prob=prob, tcc=tcc))
  if (TCC==FALSE)  return(prob)
 }

