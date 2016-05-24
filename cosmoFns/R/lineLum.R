lineLum <-
# line luminosity L in solar luminosities
# From SDR (1992), ApJ 398, L29, equation (1)
function(intInt, z, f.rest=115.27,
         omega.m = 0.27, omega.lambda = 0.73, H.0 = 71) {
    1.04e-3*intInt*f.rest/(1+z)*
        D.L(z, omega.m=omega.m, omega.lambda=omega.lambda, H.0=H.0)^2
}

