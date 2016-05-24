ordprob.pair <-
function(param,K,x,ss,data,same.means=FALSE){
      
if(!is.null(x)){
x = as.matrix(x)
p <- dim(x)[2]
}
else p = 0

      data = as.matrix(data)
      n = nrow(data)
      q = ncol(data)
      
      delta = param[1:(K-2)]
      a = delta2a(delta)
      beta =  param[(1:(p))+(K-2)]
      if(same.means){
      xi =  rep(param[K-2+p+1], q)
      corr = matrix(1,q,q)
      corr[upper.tri(corr)] = param[-c(1:(K-2+p+1))] 
      corr[lower.tri(corr)]=t(corr)[lower.tri(t(corr))]
      }
      else{
      xi =  param[(1:(q))+(K-2+p)]
      corr = matrix(1,q,q)
      corr[upper.tri(corr)] = param[-c(1:(K-2+p+q))] 
      corr[lower.tri(corr)]=t(corr)[lower.tri(t(corr))]
      }
if(!is.null(x)){
mu = x %*% matrix(beta,p,1)
}
else mu = rep(0, n)
      
      res = 0.0
      
      out = .C("cat_pair_llik_real2", 
               res = as.double(res),
               data = as.double(data), 
               corr = as.double(corr), 
               mu = as.double(mu), 
               xi = as.double(xi), 
               ss = as.double(ss),
               tresh = as.double(a),
               n = as.integer(n),
               q = as.integer(q),
               nlev = as.integer(K),
               two = as.integer(2), PACKAGE= "PLordprob")
      -out$res
   }
