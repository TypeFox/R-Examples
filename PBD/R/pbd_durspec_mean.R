pbd_durspec_mean = function(pars)
{
la3 = pars[1]
la2 = pars[2]
mu2 = pars[3]
if(mu2 == 0)
{
    rho_mean = 1/la3 * log(1 + la3/la2)
} else {
    if(la2 < Inf)
    {
       D = sqrt((la2 + la3)^2 + 2*(la2 - la3) * mu2 + mu2^2)
       rho_mean = 2/(D - la2 + la3 - mu2) * log(2 / (1 + (la2 - la3 + mu2)/D))
    } else {
       rho_mean = 0
    }
}
return(rho_mean)
}
