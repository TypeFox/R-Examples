lands.cond.t.prop.density.polar <-
function (theta, nu, zeta) 
{
    exp((nu - 1) * log(cos(theta)) - (nu/2) * log(nu) + (1 + 
        nu) * zeta * sin(theta))
}
