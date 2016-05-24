dimmingFactor <-
function(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71) {
    (1+z)/D.L(z, omega.m=omega.m, omega.lambda=omega.lambda, H.0=H.0)^2
}

