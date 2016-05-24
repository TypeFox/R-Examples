"sde.sim.shoji" <-
function(X0, t0, Dt, N, M, d1, d1.x, d1.xx, d1.t, s1){
return( .Call("sde_sim_shoji",  X0, t0, Dt, as.integer(N), as.integer(M),
 d1, d1.x, d1.xx, d1.t, s1(1,1), .GlobalEnv, PACKAGE="sde") ) 
}
