evaluate.local.maximums.lnorm <- function(x, epsilon, p, eta, sq.var, theta.fix, theta.var, x.lb, x.rb)
{
    local.extremums <- uniroot.all(function(x){num.deriv(Psi.lnorm, x, epsilon, p, eta, sq.var, theta.fix, theta.var)}, c(x.lb, x.rb))
    local.maximums <- local.extremums[second.order.num.deriv(Psi.lnorm, local.extremums, epsilon, p, eta, sq.var, theta.fix, theta.var) < 0]
    if(Psi.lnorm(x.lb, p, eta, sq.var, theta.fix, theta.var) > Psi.lnorm(x.lb + epsilon, p, eta, sq.var, theta.fix, theta.var))
        local.maximums <- c(local.maximums, x.lb)
    if(Psi.lnorm(x.rb - epsilon, p, eta, sq.var, theta.fix, theta.var) < Psi.lnorm(x.rb, p, eta, sq.var, theta.fix, theta.var))
        local.maximums <- c(local.maximums, x.rb)
    local.maximums
}