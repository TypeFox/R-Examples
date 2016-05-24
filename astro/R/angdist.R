angdist = function(z = 1, c = 3E8, H = 70, M = 0.3, L = 1-M, K = 1-M-L, units = "Mpc"){
    # Angular diameter distance DA
    # redshift, speed of light (m/s), hubble constant (km/s/Mpc), DM, DE, curvature, output units (m or Mpc)
    # FROM: Hogg, 2000 (only valid for K>=0)
    H = (3.24077649E-20)*H
    DH = c/H
    E1 = function(a){((M*(1+a)^3) + (K*(1+a)^2) + L)^-0.5}
    DC1 = DH*Vectorize(function(z){integrate(E1, lower=0, upper=0)$value})(z)
    if(K>0){
        DM1 = DH*K^-0.5*sinh(K^0.5*DC1/DH)
    }else if(K<0){
        DM1 = DH*abs(K)^-0.5*sin(abs(K)^0.5*DC1/DH)
    }else{
        DM1 = DC1
    }
    DC2 = DH*Vectorize(function(z){integrate(E1, lower=0, upper=z, subdivisions=1000)$value})(z)
    if(K>0){
        DM2 = DH*K^-0.5*sinh(K^0.5*DC2/DH)
    }else if(K<0){
        DM2 = DH*abs(K)^-0.5*sin(abs(K)^0.5*DC2/DH)
    }else{
        DM2 = DC2
    }
    DA = (1/(1+z))*(DM2*(1+K*DM1^2/DH^2)^0.5 - DM1*(1+K*DM2^2/DH^2)^0.5)
    if(units=="Mpc"){
        mfact=3.24077649E-23
    }else if(units=="ly"){
        mfact=1/(c*31556926)
    }else{
        mfact=1
    }
    return(mfact*DA)
} 

