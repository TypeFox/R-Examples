comovvol = function(z = 1, c = 3E8, H = 70, M = 0.3, L = 1-M, K = 1-M-L, units = "Gpc3"){
    # Comoving volume (all sky)
    # redshift, speed of light (m/s), hubble constant (km/s/Mpc), DM, DE, curvature, output units (m^3 or Gpc^3)
    # FROM: Hogg, 2000
    H = (3.24077649E-20)*H
    DH = c/H
    E1 = function(a){((M*(1+a)^3) + (K*(1+a)^2) + L)^-0.5}
    DC = DH*Vectorize(function(z){integrate(E1, lower=0, upper=z)$value})(z)
    if(K>0){
        DM = DH*K^-0.5*sinh(K^0.5*DC/DH)
    }else if(K<0){
        DM = DH*abs(K)^-0.5*sin(abs(K)^0.5*DC/DH)
    }else{
        DM = DC
    }
    if(K>0){
        VC = (4*pi*DH^3/(2*K))*(DM/DH*((1+K*DM^2/DH^2)^0.5)-(abs(K)^-0.5)*asinh((abs(K)^0.5)*DM/DH))
    }else if(K<0){
        VC = (4*pi*DH^3/(2*K))*(DM/DH*((1+K*DM^2/DH^2)^0.5)-(abs(K)^-0.5)*asin((abs(K)^0.5)*DM/DH))
    }else{
        VC = (4*pi/3)*(DM^3)
    }
    # units are currently in m3
    if(units=="Gpc3"){
        mfact=((3.24077649E-23)^3)/1000000000
    }else if(units=="Mpc3"){
        mfact=(3.24077649E-23)^3
    }else if(units=="ly3"){
        mfact=((1/(c*31556926))^3)/1000000000
    }else{
        mfact=1
    }
    return(mfact*VC)
}

