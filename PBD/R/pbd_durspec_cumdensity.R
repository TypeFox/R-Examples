pbd_durspec_cumdensity = function(pars,tau)
{
   integrate(function(x) pbd_durspec_density(pars,x),lower = 0,upper = tau, abs.tol = 1e-10)$value 
}