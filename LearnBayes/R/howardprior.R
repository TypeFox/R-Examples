howardprior=function (xy, par) 
{
    alpha = par[1]
    beta = par[2]
    gam = par[3]
    delta = par[4]
    sigma = par[5]
    p1 = xy[1]
    p2 = xy[2]
    u = log(p1/(1 - p1) * (1 - p2)/p2)/sigma
    z = -0.5 * u^2 + (alpha - 1) * log(p1) + (beta - 1) * log(1 - 
        p1) + (gam - 1) * log(p2) + (delta - 1) * log(1 - p2)
    return(z)
}
