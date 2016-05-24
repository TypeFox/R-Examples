"sde.sim.milstein2" <-
function(X0,  t0, Dt, N, M,
          d1, d1.x, d1.xx, s1, s1.x, s1.xx){
   return( .Call("sde_sim_milstein2",  X0,  t0, Dt, as.integer(N), 
   as.integer(M), d1, d1.x, d1.xx, s1, s1.x, s1.xx, .GlobalEnv, PACKAGE="sde") )
}

