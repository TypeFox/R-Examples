pbd_durspec_moment = function(pars,order)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
if(la2 < Inf)
{
   rho_moment = integrate(function(x) {x^order * pbd_durspec_density(pars,x) },lower = 0, upper = Inf, abs.tol = 1e-10)
} else {
   if(order == 0)
   {
       rho_moment = 1
   } else {
       rho_moment = 0
   }
}
return(rho_moment$value)
}
