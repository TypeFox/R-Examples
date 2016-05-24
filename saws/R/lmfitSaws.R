`lmfitSaws`<-function(x,y){
    K<-length(y)
    p<-dim(x)[[2]]
    beta<- solve(t(x) %*% x ) %*% t(x) %*% matrix(y,K,1)
    resid<- y - x %*% beta
    omega<-array(NA,dim=c(K,p,p))
    u<-matrix(NA,K,p)
    for (i in 1:K){
        omega[i,,]<- matrix(x[i,],p,1) %*% matrix(x[i,],1,p) 
        u[i,]<- x[i,]*resid[i]
    }
    out<-list(coefficients=beta,u=u,omega=omega)
    out
}
