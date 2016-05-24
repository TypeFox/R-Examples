Ubeta.fun <-
function(ord.delta,ord.z,ord.bz,Rt,ord.wt,r){
    if(r == 0) a = ord.z*ord.wt*(ord.delta-Rt*exp(-ord.bz))
    if(r > 0) a = ord.z*ord.wt*(ord.delta-(1/r)*log(1+r*Rt*exp(-ord.bz)))
    b = apply(a,2,sum)
    return(b)
}
