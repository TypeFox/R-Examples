hydraulicCutOff <- 
function(theta_R, k0, k1, n, x0 = 6.653)
{
    exp((1/n) * log(-k0 / (log(exp(-k0/x0^n) - theta_R/k1))))
}
