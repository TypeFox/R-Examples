Hbeta.fun <-
function(ord.z,ord.bz,Rt,ord.wt,r,p){
    a = matrix(0,nrow=p,ncol=p)
    for(j in 1:p){
      for(k in 1:p){
         a[j,k] = -sum(ord.wt*as.numeric(ord.z[,j])*as.numeric(ord.z[,k])*Rt*exp(-ord.bz)/(1+r*Rt*exp(-ord.bz)))
      }
    }
    return(a)
}
