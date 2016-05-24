age = function(z = 1, H = 70, M = 0.3, L = 1-M, K = 1-M-L, units = "Gyr"){
    # Age of the universe
    # redshift, hubble constant (km/s/Mpc), DM, DE, curvature, output units (s or Gyr)
    # FROM: Hogg, 2000
    H = (3.24077649E-20)*H
    tH = 1/H
    E2 = function(a){ ((1+a)^-1.0)*(((M*(1+a)^3) + (K*(1+a)^2) + L)^-0.5) }
    tL = tH*Vectorize(function(z){integrate(E2, lower=z, upper=Inf)$value})(z)
    if(units=="Gyr"){
        mfact=3.16887646E-17
    }else{
        mfact=1
    }
    return(mfact*tL)
}

