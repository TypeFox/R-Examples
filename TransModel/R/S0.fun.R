S0.fun <-
function(ord.bz,Rt,ord.wt,r,n){
    Rt0 = c(0,Rt[1:(n-1)])
    s0 = numeric(n)
    if(r == 0){
       a = ord.wt*exp(-ord.bz)
       s0 = sum(a) - cumsum(a) + a
    }
    if(r > 0){
       for(i in 1:n){
          a = ord.wt*exp(-ord.bz)/(1+r*Rt0[i]*exp(-ord.bz))
          s0[i] = sum(a[i:n])
       }
    }
    return(s0)
}
