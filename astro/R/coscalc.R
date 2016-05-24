coscalc = function(z = 1, c = 3E8, H = 70, M = 0.3, L = 1-M, K = 1-M-L, dunit = "Mpc", vunit = "Gpc3", tunit = "Gyr", r = 1, inp = "arcsec", out = "kpc"){
    
    # fake top-level funciton
    temp = function(z, c, H, M, L, K, dunit, vunit, tunit, r, inp, out){
    
    # units
    if(dunit=="Mpc"){
        dfact=3.24077649E-23
    }else if(dunit=="ly"){
        dfact=1/(c*31556926)
    }else{
        dfact=1
    }
    
    if(vunit=="Gpc3"){
        vfact=((3.24077649E-23)^3)/1000000000
    }else if(vunit=="ly3"){
        vfact=((1/(c*31556926))^3)/1000000000
    }else{
        vfact=1
    }
    
    if(tunit=="Gyr"){
        tfact=3.16887646E-17
    }else{
        tfact=1
    }
    
    # calculations
    H = (3.24077649E-20)*H
    DH = c/H
    E1 = function(a){((M*(1+a)^3) + (K*(1+a)^2) + L)^-0.5}
    DC = DH*Vectorize(function(z){integrate(E1, lower=0, upper=z, subdivisions=1000)$value})(z)
    codist.los = dfact*DC       # codist.los
    if(K>0){
        DM = DH*K^-0.5*sinh(K^0.5*DC/DH)
        VC = (4*pi*DH^3/(2*K))*(DM/DH*((1+K*DM^2/DH^2)^0.5)-(abs(K)^-0.5)*asinh((abs(K)^0.5)*DM/DH))
    }else if(K<0){
        DM = DH*abs(K)^-0.5*sin(abs(K)^0.5*DC/DH)
        VC = (4*pi*DH^3/(2*K))*(DM/DH*((1+K*DM^2/DH^2)^0.5)-(abs(K)^-0.5)*asin((abs(K)^0.5)*DM/DH))
    }else{
        DM = DC
        VC = (4*pi/3)*(DM^3)
    }
    codist.trans = dfact*DM     # codist.trans
    DM1 = 0
    DM2 = DM
    DA = (1/(1+z))*(DM2*(1+K*DM1^2/DH^2)^0.5 - DM1*(1+K*DM2^2/DH^2)^0.5)
    angdist = dfact*DA          # angdist
    DL = (1+z)*DM
    lumdist = dfact*DL          # lumdist
    covol = vfact*VC            # covol
    tH = 1/H
    E2 = function(a){ ((1+a)^-1.0)*(((M*(1+a)^3) + (K*(1+a)^2) + L)^-0.5) }
    tL1 = tH*Vectorize(function(z){integrate(E2, lower=0, upper=z)$value})(z)
    lookback = tfact*tL1        # lookback
    tL2 = tH*Vectorize(function(z){integrate(E2, lower=z, upper=Inf)$value})(z)
    age = tfact*tL2             # age
    
    # angular size
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
    angsize = result            # angsize
    
    # return results
    return=c(z=z, codist.los=codist.los, codist.trans=codist.trans, angdist=angdist, lumdist=lumdist, covol=covol, lookback=lookback, age=age, angsize=angsize)
    
    }

    return(t(Vectorize(temp)(z, c, H, M, L, K, dunit, vunit, tunit, r, inp, out)))

}

