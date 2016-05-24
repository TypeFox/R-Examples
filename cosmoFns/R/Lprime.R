Lprime <-
# Line luminosity L' in K km/s pc^-2
# From SDR (1992), ApJ 398, L29, equation (3)
function(intInt, z, f.rest=115.27, omega.m = 0.27, omega.lambda = 0.73,
         H.0 = 71) {
    3.25e7*intInt*
        D.L(z, omega.m=omega.m, omega.lambda=omega.lambda, H.0=H.0)^2/
            (f.rest^2*(1+z))
}

