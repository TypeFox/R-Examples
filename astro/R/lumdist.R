lumdist = function(z = 1, c = 3E8, H = 70, M = 0.3, L = 1-M, K = 1-M-L, units = "Mpc"){
    # Luminosty distance DL
    # redshift, speed of light (m/s), hubble constant (km/s/Mpc), DM, DE, curvature, output units (m or Mpc)
    # FROM: Hogg, 1999
    H = (3.24077649E-20)*H
    DH = c/H
    E1 = function(a){((M*(1+a)^3) + (K*(1+a)^2) + L)^-0.5}
    DC = DH*Vectorize(function(z){integrate(E1, lower=0, upper=z)$value})(z)
    if(K>0){
        DM = DH*(K^-0.5)*sinh(K^0.5*DC/DH)
    }else if(K<0){
        DM = DH*(abs(K)^-0.5)*sin((abs(K)^0.5)*DC/DH)
    }else{
        DM = DC
    }
    DL = (1+z)*DM
    if(units=="Mpc"){
        mfact=3.24077649E-23
    }else if(units=="ly"){
        mfact=1/(c*31556926)
    }else{
        mfact=1
    }
    return(mfact*DL)
}

