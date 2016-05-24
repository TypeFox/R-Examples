"sde.sim.KPS" <-
function(X0,  t0, Dt, N, M, d1, d1.x, d1.xx, s1, s1.x, s1.xx, Z, U){
    return( .Call("sde_sim_KPS",  X0,  t0, Dt, as.integer(N), as.integer(M),
	d1, d1.x, d1.xx, s1, s1.x, s1.xx, Z, U, .GlobalEnv, PACKAGE="sde") )
}

