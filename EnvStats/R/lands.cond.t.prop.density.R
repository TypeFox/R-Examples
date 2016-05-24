lands.cond.t.prop.density <-
function (tau, nu, zeta) 
{
    ((nu + tau^2)^(-(nu + 1)/2)) * exp(((nu + 1) * zeta * tau)/sqrt(nu + 
        tau^2))
}
