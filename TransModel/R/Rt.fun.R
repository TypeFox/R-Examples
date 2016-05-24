Rt.fun <-
function(ord.delta,ord.wt,s0){
    Rt = cumsum(ord.wt*ord.delta/s0)
    return(Rt)
}
