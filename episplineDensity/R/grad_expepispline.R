grad_expepispline <- function (epiparameters, r, x, moment=0) 
{
#
# Gradient for exponential epispline.
#
n <- length(x)
e.grad <- matrix(0, length (r), n)
for (i in 1:n) {
    if (x[i] <= epiparameters$mN && x[i] >= epiparameters$m0)
    {
        c.out <- coeff (x[i], epiparameters)
        e.grad[,i] <- exp (sum (- c.out * r)) * c.out # positive for now
        e.grad[,i] <- e.grad[,i] * (-1) * (x[i] ^ moment)
    }
}
return (e.grad)
}
