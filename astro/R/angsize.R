angsize = function(z = 1, r = 1, inp = "arcsec", out = "kpc", c = 3E8, H = 70, M = 0.3, L = 1-M, K = 1-M-L){
    # Apparent Angular Size RA
    # redshift, radius of object, speed of light (m/s), hubble constant (km/s/Mpc), DM, DE, curvature, input units (m, pc, kpc or Mpc), output units (deg, rad, arcsec)
    # FROM: Hogg, 2000 (only valid for K>=0)
    inverse=FALSE
    if(inp=="deg"|inp=="rad"|inp=="arcsec"){
        inverse=TRUE
        real_r=r
        real_inp=inp
        real_out=out
        r=1
        inp="kpc"
        out=real_inp
    }
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
    if(inp=="pc"){
        mfact2=3.24077649E-17
    }else if(inp=="kpc"){
        mfact2=3.24077649E-20
    }else if(inp=="Mpc"){
        mfact2=3.24077649E-23
    }else{
        mfact2=1
    }
    DA=mfact2*DA
    theta=atan(r/DA) # radians
    if(out=="arcsec"){
        mfact=206264.806
    }else if(out=="deg"){
        mfact=57.2957795
    }else{
        mfact=1
    }
    result=theta*mfact
    if(inverse){
        result=real_r*(1/result)
        if(real_out=="m"){
            mfact=3.08568025E+19
        }else if(real_out=="pc"){
            mfact=10^3
        }else if(real_out=="Mpc"){
            mfact=10^-3
        }else{
            mfact=1
        }
        result=result*mfact
    }
    return(result)
}

