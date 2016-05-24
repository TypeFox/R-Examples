evaluate.local.maximums <- function(x, epsilon, p, eta, theta.fix, theta.var, x.lb, x.rb)
{
    local.extremums <- uniroot.all(function(x){num.deriv(Psi, x, epsilon, p, eta, theta.fix, theta.var)}, c(x.lb, x.rb))
    local.maximums <- local.extremums[second.order.num.deriv(Psi, local.extremums, epsilon, p, eta, theta.fix, theta.var) < 0]
    if(Psi(x.lb, p, eta, theta.fix, theta.var) > Psi(x.lb + epsilon, p, eta, theta.fix, theta.var))
        local.maximums <- c(local.maximums, x.lb)
    if(Psi(x.rb - epsilon, p, eta, theta.fix, theta.var) < Psi(x.rb, p, eta, theta.fix, theta.var))
        local.maximums <- c(local.maximums, x.rb)
    local.maximums
}

