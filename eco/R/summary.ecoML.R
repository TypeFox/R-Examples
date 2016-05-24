##for simlicity, this summary function only reports parameters related to W_1 and W_2
summary.ecoML <- function(object, CI = c(2.5, 97.5),  param = TRUE, units = FALSE, subset = NULL, ...) { 


      n.col<-5
      if(object$context) n.col<-7
      if (object$fix.rho) n.col<-n.col-1


  n.row<-1
  if (object$sem) n.row<-3


  param.table<-matrix(NA, n.row, n.col)
   
  if (!object$context) {
   param.table[1,]<-object$theta.em 
   cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho")
  }
  else if (object$context && !object$fix.rho) {
   cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho1X","rho2X","rho12")
   param.table[1,]<-object$theta.em[c(2,3,5,6,7,8,9)]   
  }
  else if (object$context && object$fix.rho) {
   cname<-c("mu1", "mu2", "sigma1", "sigma2", "rho1X","rho2X")
   param.table[1,]<-object$theta.em[c(2,3,5,6,7,8)] 
  }

  if (n.row>1) {
    if (!object$context) {
    param.table[2,]<-sqrt(diag(object$Vobs))
    param.table[3,]<-Fmis<-1-diag(object$Iobs)/diag(object$Icom)
   }
   else if (object$context && !object$fix.rho) {
    param.table[2,]<-sqrt(diag(object$Vobs)[c(2,3,5,6,7,8,9)])
    param.table[3,]<-Fmis<-(1-diag(object$Iobs)/diag(object$Icom))[c(2,3,5,6,7,8,9)]
  }
   else if (object$context && object$fix.rho) {
    param.table[2,]<-sqrt(diag(object$Vobs)[c(2,3,5,6)])
    param.table[3,]<-Fmis<-(1-diag(object$Iobs)/diag(object$Icom))[c(2,3,5,6)]
  }

  }
  rname<-c("ML est.", "std. err.", "frac. missing")
  


  rownames(param.table)<-rname[1:n.row]
  colnames(param.table)<-cname[1:n.col]
  
  n.obs <- nrow(object$W)

   if (is.null(subset)) subset <- 1:n.obs 
  else if (!is.numeric(subset))
    stop("Subset should be a numeric vector.")
  else if (!all(subset %in% c(1:n.obs)))
    stop("Subset should be any numbers in 1:obs.")
    
  table.names<-c("mean", "std.dev", paste(min(CI), "%", sep=" "),
                 paste(max(CI), "%", sep=" ")) 

  W1.mean <- mean(object$W[,1])
  W2.mean <- mean(object$W[,2])
  W1.sd <- sd(object$W[,1])
  W2.sd <- sd(object$W[,2])
#  W1.q1 <-  W1.mean-1.96*W1.sd
#  W1.q2 <-  W1.mean+1.96*W1.sd
#  W2.q1 <-  W2.mean-1.96*W2.sd
#  W2.q2 <-  W2.mean+1.96*W2.sd
  W1.q1 <-  quantile(object$W[,1],min(CI)/100)
  W1.q2 <-  quantile(object$W[,1],max(CI)/100)
  W2.q1 <-  quantile(object$W[,2],min(CI)/100)
  W2.q2 <-  quantile(object$W[,2],max(CI)/100)
  
  agg.table <- rbind(cbind(W1.mean, W1.sd, W1.q1, W1.q2),
                     cbind(W2.mean, W2.sd, W2.q1, W2.q2)) 
  colnames(agg.table) <- table.names
  rownames(agg.table) <- c("W1", "W2")

 # if (is.null(object$N))
 #   N <- rep(1, nrow(object$X))
 # else
 
 agg.wtable<-NULL
 if (!is.null(object$N)) {
    N <- object$N
}
else {
    N <- rep(1:n.obs)
}
  weighted.var <- function(x, w) {
    return(sum(w * (x - weighted.mean(x,w))^2)/((length(x)-1)*mean(w)))
    }

  W1.mean <- weighted.mean(object$W[,1], object$X*N)
  W2.mean <- weighted.mean(object$W[,2], (1-object$X)*N)
  W1.sd <- weighted.var(object$W[,1], object$X*N)^0.5
  W2.sd <- weighted.var(object$W[,1], (1-object$X)*N)^0.5
  W1.q1 <-  W1.mean-1.96*W1.sd
  W1.q2 <-  W1.mean+1.96*W1.sd
  W2.q1 <-  W2.mean-1.96*W2.sd
  W2.q2 <-  W2.mean+1.96*W2.sd
#  W1.q1 <-  quantile(object$W[,1] * object$X*N/mean(object$X*N),min(CI)/100)
#  W1.q2 <-  quantile(object$W[,1] * object$X*N/mean(object$X*N),max(CI)/100)
#  W2.q1 <-  quantile(object$W[,2]*(1-object$X)*N/(mean((1-object$X)*N)),min(CI)/100)
#  W2.q2 <-  quantile(object$W[,2]*(1-object$X)*N/(mean((1-object$X)*N)),max(CI)/100)
  agg.wtable <- rbind(cbind(W1.mean, W1.sd, W1.q1, W1.q2),
                      cbind(W2.mean, W2.sd, W2.q1, W2.q2))
  colnames(agg.wtable) <- table.names
  rownames(agg.wtable) <- c("W1", "W2")
  
  
  if (units) 
    W.table <- object$W[subset,] 
  else
    W.table <-  NULL
  
  ans <- list(call = object$call, iters.sem = object$iters.sem,
              iters.em = object$iters.em, epsilon = object$epsilon,
              sem = object$sem, fix.rho = object$fix.rho, loglik = object$loglik,
              rho=object$rho, param.table = param.table, W.table = W.table, 
              agg.wtable = agg.wtable, agg.table=agg.table, n.obs = n.obs) 
 # if (object$fix.rho)
 #   ans$rho<-object$rho
  
  class(ans) <-"summary.ecoML"
  return(ans)
}
