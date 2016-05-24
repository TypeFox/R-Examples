"sde.sim.milstein" <-
function(X0,  t0, Dt, N, M, d1, s1, s1.x){
   return( .Call("sde_sim_milstein",  X0,  t0, Dt, as.integer(N), 
   as.integer(M), d1, s1, s1.x, .GlobalEnv, PACKAGE="sde") )
}

