extract_hica <-
function(energy.obj,comp,level){
 components=energy.obj$components
 energy=energy.obj$energy
 HICA.obj=energy.obj$HICA.obj
 X=HICA.obj$X
 bas=HICA.obj$basis
 maxlev=dim(HICA.obj$aggregation)[1]
 p <- dim(X)[2]
 if (comp<1 || comp > p){
  stop(paste("the extracted component must to be between 1 and ",p,sep=""))
 }
 if (level<1 || level > maxlev){
  stop(paste("the level must to be between 1 and ",maxlev,sep=""))
 }
 loading=bas[[level]][,components[1:comp,level]]
 S=X%*%loading%*%solve(t(loading)%*%loading)
 cum.energy=energy[,level]
 list(X=X,S=S,C=loading, cum.energy=cum.energy)
}
