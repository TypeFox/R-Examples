mass.CO <-
function(intInt, z, alpha=0.8, f.rest=115.27, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71) {
    Lprime(intInt, z, f.rest,
           omega.m=omega.m, omega.lambda=omega.lambda, H.0=H.0) * alpha
}

