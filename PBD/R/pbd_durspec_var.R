pbd_durspec_var = function(pars)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
rho_var = pbd_durspec_moment(pars,2) - (pbd_durspec_mean(pars))^2
return(rho_var)
}
