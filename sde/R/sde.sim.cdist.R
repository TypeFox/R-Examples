"sde.sim.cdist" <-
function(X0, t0, Dt, N, M, rcdist=NULL, theta=NULL){
 cds <-function(t,x) rcdist(1,t,x,theta)
 return( .Call("sde_sim_cdist",  X0, t0, Dt, as.integer(N), as.integer(M),
  cds, .GlobalEnv, PACKAGE="sde") ) 
}

