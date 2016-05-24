CutoffW <-
function(res,delta,cu=1.8554,zmax=10,gridsize=2000){
# computes cutoff on residual scale 
xs <- seq(from=cu,to=zmax,length=gridsize)
fs <- apply(as.matrix(xs),1,RappW,res=res,delta=delta)
alpha <- min(min(fs[2,]),1)
xu <- max(xs[fs[1,] < alpha])
xl <- Izero(xu)
list(alpha=alpha,tl=xl,tu=xu)}

