"sde.sim.euler" <-
function(X0, t0, Dt, N, M, d1, s1, s1.x, alpha, eta, pred.corr){
 return( .Call("sde_sim_euler",  X0, t0, Dt, as.integer(N), as.integer(M),
 d1, s1, s1.x, alpha, eta, as.logical(pred.corr), .GlobalEnv, PACKAGE="sde") ) 
}

