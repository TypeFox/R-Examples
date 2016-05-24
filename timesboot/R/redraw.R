redraw <-
function(series,...){
args <- list(...)
ws = spec.pgram(series,spans=args$p,plot=FALSE,detrend=args$detrend)
return(ws$spec)
}
