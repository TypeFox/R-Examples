schechter = function(x, knee, slope, norm, bw = 0.1, mag = FALSE, log = FALSE, ...){
    
    # driver 2008 (magnitude)
    # bw = 0.5, mag = TRUE, log = FALSE, knee = -21.32, slope = -1.32, norm = 4.8E-3
    # mag = seq(-24,-17,len=100)
    # aplot(mag, 0.5*schechter(mag, bw = 0.5, mag = TRUE, log = FALSE, knee = -21.32, slope = -1.32, norm = 4.8E-3), log="y", yformat="p", las=1, ylim=c(1e-6,1e-2), type="l", xaxs="i", yaxs="i", nxmin=1)
    
    # baldry 2012 (stellar mass)
    # bw = 0.1, mag = FALSE, log = TRUE, knee = 10^10.66, slope = c(-0.35,-1.47), norm = c(3.96E-3,0.79E-3)
    # mass = 10^seq(7,11.6,len=100)
    # aplot(mass, schechter(mass, bw = 0.1, mag = FALSE, log = TRUE, knee = 10^10.66, slope = c(-0.35,-1.47), norm = c(3.96E-3,0.79E-3)), log="xy", yformat="p", las=1, type="l", xaxs="i", yaxs="i", ylim=c(1e-5,2e-1))
    
    # check inputs have same length (for double schechter fits)
    if(any(is.na(knee))){knee = knee[-which(is.na(knee))]}
    if(any(is.na(slope))){slope = slope[-which(is.na(slope))]}
    if(any(is.na(norm))){norm = norm[-which(is.na(norm))]}
    lengths = c(length(knee), length(slope), length(norm))
    if(max(lengths)!=min(lengths)){
        knee = rep(knee,max(lengths))[1:max(lengths)]
        slope = rep(slope,max(lengths))[1:max(lengths)]
        norm = rep(norm,max(lengths))[1:max(lengths)]
    }
    
    # mag check
    if(mag){
        log = FALSE
    }
    
    # single schechter
    funclin1 = function(x, k, s, n){
        phi = n/k * ((x/k)^s) * exp(-x/k)
        return(phi)
    }
    
    # double schechter
    funclin2 = function(x, k1, s1, n1, k2, s2, n2){
        phi = (n1/k1 * ((x/k1)^s1) * exp(-x/k1)) + (n2/k2 * ((x/k2)^s2) * exp(-x/k2))
        return(phi)
    }
    
    # single schechter (log)
    funclog1 = function(x, k, s, n){
        phi = log(10) * n * (10^((s+1)*(x-k))) * exp(-10^(x-k))
        return(phi)
    }
    
    # double schechter (log)
    funclog2 = function(x, k1, s1, n1, k2, s2, n2){
        phi = (log(10) * n1 * (10^((s1+1)*(x-k1))) * exp(-10^(x-k1))) + (log(10) * n2 * (10^((s2+1)*(x-k2))) * exp(-10^(x-k2)))
        return(phi)
    }
    
    # single schechter (magnitude)
    funcmag1 = function(x, k, s, n){
        phi = 0.4 * log(10) * n * 10^(-0.4*(x-k)*(s+1)) * exp(-10^(-0.4*(x-k)))
        return(phi)
    }
    
    # double schechter (magnitude)
    funcmag2 = function(x, k1, s1, n1, k2, s2, n2){
        phi = (0.4 * log(10) * n1 * 10^(-0.4*(x-k1)*(s1+1)) * exp(-10^(-0.4*(x-k1)))) + (0.4 * log(10) * n2 * 10^(-0.4*(x-k2)*(s2+1)) * exp(-10^(-0.4*(x-k2))))
        return(phi)
    }
    
    # integration steps
    stepslo = x-bw/2
    stepshi = x+bw/2
    
    # integrate
    phi = {}
    for(i in 1:length(x)){
        
        if(mag){
            # magnitudes
            if(max(lengths)==1){
                temp = integrate(funcmag1, lower=stepslo[i], upper=stepshi[i], k=knee[1], s=slope[1], n=norm[1], ...)
            }else{
                temp = integrate(funcmag2, lower=stepslo[i], upper=stepshi[i], k1=knee[1], s1=slope[1], n1=norm[1], k2=knee[2], s2=slope[2], n2=norm[2], ...)
            }
        }else if(log){
            # logged
            if(max(lengths)==1){
                temp = integrate(funclog1, lower=stepslo[i], upper=stepshi[i], k=knee[1], s=slope[1], n=norm[1], ...)
            }else{
                temp = integrate(funclog2, lower=stepslo[i], upper=stepshi[i], k1=knee[1], s1=slope[1], n1=norm[1], k2=knee[2], s2=slope[2], n2=norm[2], ...)
            }
        }else{
            # linear
            if(max(lengths)==1){
                temp = integrate(funclin1, lower=stepslo[i], upper=stepshi[i], k=knee[1], s=slope[1], n=norm[1], ...)
            }else{
                temp = integrate(funclin2, lower=stepslo[i], upper=stepshi[i], k1=knee[1], s1=slope[1], n1=norm[1], k2=knee[2], s2=slope[2], n2=norm[2], ...)
            }
        }
        phi = c(phi, temp$value)
        
    }
    
    # return data
    return(phi/bw)
    
}

