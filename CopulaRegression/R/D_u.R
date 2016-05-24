D_u <-
function(u,v,theta,family){
    u[u<0]=0
    u[u>1]=1
    v[v<0]=0
    v[v>1]=1
   out<-BiCopHfunc(u,v,family=family,par=theta)$hfunc1
   return(out)
}
