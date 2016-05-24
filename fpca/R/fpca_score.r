#########
fpca.score<-function(data.m,grids.u,muhat,eigenvals,eigenfuncs,sig2hat,K){
##estimated conditional principal component scores (BLUPs): \hat{E(\xi_k|Y)}, k=1,...,K
##Name:FPcScore
##para:
##     data.m -- data matrix; same as input for fpca.mle
##     grids.u -- grid of time points used in evaluating the mean and eigenfunctions (on the original scale); (returned by fpca. mle)
##     muhat,eigenvals, eigenfuncs, sig2hat -- (estimated) mean, eigenvalues, eigenfunctions and noise variance; (returned by fpca.mle)
##     K -- number of eigenfunctions used in the model, i.e., (estimated) dimension of the process
##return: first K conditional PC scores (the BLUP estimates): n by K
temp<-table(data.m[,1])
n<-length(temp)             ##     number of curves;
m.l<-as.vector(temp)        ##     m.l -- number of time points per curve
result<-matrix(0,n,K)       ##First K FPC scores for each subject

N <- length(grids.u)        ## number of time points on the grid
evalmat <- diag(eigenvals[1:K])  ## diagonal matrix of the first K (estimated) eigenvalues
current <- 0  ## current index
eigenfuncs.u<-t(eigenfuncs)   ## dimmension: grid_length by K

data.u<-matrix(as.numeric(as.vector(data.m[,-1])),nrow=nrow(data.m[,-1]),ncol=ncol(data.m[,-1]))     ##convert obs matrix to be numierc

  for (i in 1:n){
      Y <- as.vector(data.u[(current+1):(current+m.l[i]),1])  ## observed  measurements of ith curve
      meastime <- data.u[(current+1):(current+m.l[i]),2] ## measurement times of the ith curve
      gridtime <- ceiling(N*meastime)   ## project measurement time onto the grid
      muy <- muhat[gridtime]
      Phiy  <- matrix(eigenfuncs.u[gridtime,1:K],ncol=K)
      Sigy <- Phiy %*% evalmat %*% t(Phiy) + sig2hat * diag(m.l[i])
      temp.y<-matrix(Y-muy)
      result[i,] <- evalmat %*% t(Phiy) %*% solve(Sigy,temp.y)
      current <- current + m.l[i]
  }
return(result)
}