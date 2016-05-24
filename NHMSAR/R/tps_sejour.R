tps_sejour <-
function(x,s=.0){
	w = which(x<s)
	dw = c(2,diff(w),2)
    wdw=which(dw>1)
    dur = diff(wdw)
    return(dur)	
}
