epispline <- function (epiparameters, r, x, exponentiate=FALSE, moment)
{
# epiparameters: list of epiparameters
#             r: numeric vector of coefficients
#             x: numeric vector of data
#  exponentiate: logical: compute exponential epispline?
#        moment: numeric, 1 or 2 for first or second moment. If missing
#                use 0. Moment requires exponential to be TRUE.
#
n <-  length(x) 
s <- numeric (n)
#
# For "moment" computations, multiply by x^m. If no moment is asked
# for, use x^0 for convenience.
#

if (missing(moment))
    moment <- 0
else {
    if (exponentiate == FALSE)
        stop ("Don't supply moment without exponentiate in epispline")
    if (length (moment) != 1 || !is.numeric (moment) || 
        !is.element (moment, c(1, 2)))
        stop ("Illegal moment passed to epispline; should be 1 or 2")
}
#
# For each point, compute the vector of coefficients for that point, c,
# and the compute <c . r>. For exponentiate, compute exp (- that);
# for moments, multiply by x[i]^moment.
#
for (i in 1:n) {
    if (x[i] <= epiparameters$mN && x[i] >= epiparameters$m0) {
        c.out <- coeff(x[i], epiparameters) 
        if (exponentiate)
            s[i] <- (x[i]^moment) * exp (- sum (c.out * r))
        else 
            s[i] <- sum (c.out * r)
    }
    else 
        s[i] <- Inf
}

return (s)

}
