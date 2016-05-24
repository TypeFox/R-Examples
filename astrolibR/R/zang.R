#source('lumdist.R')

zang = function(dl,z, h0, k, lambda0, omega_m, q0) {

    d = lumdist(z, h0 = h0, k = k, lambda0 = lambda0, omega_m = omega_m,  q0 = q0)

    radeg = 180/pi
    return(radeg*3600*dl*(1.+z)^2/(1000.*d))
}
