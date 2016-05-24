ordprob.mple <-
function(y, x=NULL, K, start.par=list(type="default"), same.means=FALSE, eval.max=1000, iter.max=600, ...){
      n <- nrow(y)
      q <- ncol(y)
      q.mean <- ifelse(same.means,1,q)
      if(!is.null(x)){ p <- dim(x)[2]}
      else{ p = 0}
      y <- data.matrix(y)
      vy <- var(c(y))

      if(start.par$type=="default")
      {
      start_cor = cor(y)[upper.tri(cor(y))]
      if(p>0) {start_beta = rep(0,p)}
	else{ start_beta = NULL}
      start_xi = rep(0,q.mean)
      start_thresh = rep(0,(K-2))
      start <- c(start_thresh, start_beta, start_xi, start_cor)
      }
      if(start.par$type=="list")
      {
      start_cor = start.par$cor
      start_beta = start.par$beta
      start_xi = start.par$xi
      start_thresh = start.par$thresh
      start <- c(start_thresh, start_beta, start_xi, start_cor)
      }
      low <- c(rep(-Inf,K-2), rep(-Inf,p+q.mean), rep(-0.99,q*(q-1)/2))
      upp <- c(rep(Inf,K-2),  rep(Inf,p+q.mean),  rep(0.99,q*(q-1)/2))
      
      res <- nlminb(start, ordprob.pair, lower=low, upper=upp, K=K, x=x, ss=vy, data=y, same.means=same.means, control=list(iter.max=iter.max, eval.max=eval.max))
      if(!is.null(x)){ 
	cor.out <- matrix(NA, q,q)
	cor.out[upper.tri(cor.out)] <- res$par[-c(1:(K-2+p+q.mean))]
	out <- list(thresh=res$par[1:(K-2)], beta=res$par[(1:p)+(K-2)], xi = res$par[(1:(q.mean))+(K-2+p)], cor=cor.out, opt=res)
	}      
	else{ 
	cor.out <- matrix(NA, q,q)
	cor.out[upper.tri(cor.out)] <- res$par[-c(1:(K-2+q.mean))]
	out <- list(thresh=res$par[1:(K-2)], xi = res$par[(1:(q.mean))+(K-2)],cor=cor.out, opt=res)
	}
  out
}
