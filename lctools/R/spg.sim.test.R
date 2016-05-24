spg.sim.test <- function(realizations=1000, nperm = 99){

  .Random.seed
  set.seed(10)

  #nrp = np.random.permutation
  nsqrts <-c(5, 7, 9)
  rhos <- 0.1 * (0:9)
  sigmas <- c(1.0, 2.0, 3.0, 4.0)
  alpha = 0.05
  
  res <- replicate(length(nsqrts),matrix(0,length(rhos),length(sigmas)), FALSE)
  
  gini_mean <- res
  gini_std <- res
  moran_mean <- res
  moran_res <- res
  
  pivot <- (1+nperm) * alpha
  
  for (nsqrt in nsqrts){
    i<-which(nsqrts %in% nsqrt)
    n = nsqrt^2
    
    #x = np.arange(1,n+1)
    w = lat2w(nsqrt,nsqrt)[[2]]
    
    for(rho in rhos){
      ri<-which(rhos %in% rho)
      
      A <- diag(n) - rho * w
      AI <- ginv(A)
      
      for(sigma in sigmas){
        
        si<-which(sigmas %in% sigma)
        
        print (paste(nsqrt,n, sigma, rho,sep=" - "))
        
        gm <- 0
        gv <- 0
        ginis <- rep(0,realizations)
        morans <- rep(0,realizations)
        counter <- 0
        
        for(realization in 1:realizations){
          e <- rnorm(n)*sigma
          x <-AI %*% e
          x <- x + abs(min(x)) # make nonnegative
               
          res0 <- spGini.w(x, w)
          gini <- res0[[1]]
          ginis[counter] <- gini
          mi <- moransI.w(x,w)
          morans[counter] <- mi[[1]]
          counter <- counter + 1
          nni <- res0[[3]] # non-neighbor inequality
          
          ## Inner permutations
          mi_count <- 0
          nni_count <- 0
          for (perm in 1:nperm){
            xr <- sample(x)
            mir <- moransI.w(xr,w)[[1]]
            nni_r <- spGini.w(xr, w)[[3]] 
            if (mir >= mi[[1]]){
              mi_count <-mi_count + 1
            }
            if (nni_r >= nni){
              nni_count <- nni_count + 1
            }
          }
          if(nni_count <= pivot) {res[[i]][ri,si] <- res[[i]][ri,si] + 1}
          if(mi_count <= pivot) {moran_res[[i]][ri,si] <- moran_res[[i]][ri,si] + 1} 
          
#           if ((nperm-nni_count)< nni_count) {nni_count=nperm-nni_count}
#           if ((nperm-mi_count)< mi_count) {mi_count=nperm-mi_count}
#           

         # print(paste(realization,n,rho,sigma, "nni sig: ", res[[i]][ri,si], "I nsig: ", moran_res[[i]][ri,si],sep=" , "))
        }
        gini_mean[[i]][ri,si] = mean(ginis)
        moran_mean[[i]][ri,si] = mean(morans)
        gini_std[[i]][ri,si] = sd(ginis)
      }
    }
  }
  
  # res has the number of realizations where significant clustering is
  # obtained first dimension = sample size, second = rho, third = sigma
  
  res<-lapply(res,function(x) {x/realizations}) #res = res/(realizations*1.)
  moran_res<-lapply(moran_res,function(x) {x/realizations}) #res = res/(realizations*1.)
  
  return(list(res=res, moran_res=moran_res, gini_mean=gini_mean, gini_std = gini_std, moran_mean = moran_mean))
}


