"sde.sim.ozaki" <-
function(X0, t0, Dt, N, M, d1, d1.x, s1){
 return( .Call("sde_sim_ozaki",  X0, t0, Dt, as.integer(N), as.integer(M),
 d1, d1.x, s1(1,1), .GlobalEnv, PACKAGE="sde") ) 
}
