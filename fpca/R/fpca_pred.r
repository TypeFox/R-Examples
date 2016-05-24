
fpca.pred<-function(fpcs, muhat,eigenfuncs){
## predicted trajectory for each sample curve
##para: fpcs: (estimated) FPC score;     (returned by fpca.score)
##      muhat, eigenfuncs: (estimated) mean and eigenfunctions evaluated on a fine grid.    (returned by fpca.mle)
##return: predicted trajectories: grid_length by n
eigenfuncs.u<-t(eigenfuncs)   ## dimmension: grid_length by K
result<-muhat+eigenfuncs.u%*%t(fpcs)
return(result)                      ##each column corresponds to a predicted trajectory
}
