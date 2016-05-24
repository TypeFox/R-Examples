gauss = function(x, mean = 0, sigma = 1, norm = 1, bw = 0.1, ftype = "lin", ...){
    
    # check inputs have same length (for double gauss fits)
    if(any(is.na(mean))){mean = mean[-which(is.na(mean))]}
    if(any(is.na(sigma))){sigma = sigma[-which(is.na(sigma))]}
    if(any(is.na(norm))){norm = norm[-which(is.na(norm))]}
    lengths = c(length(mean), length(sigma), length(norm))
    if(max(lengths)!=min(lengths)){
        mean = rep(mean,max(lengths))[1:max(lengths)]
        sigma = rep(sigma,max(lengths))[1:max(lengths)]
        norm = rep(norm,max(lengths))[1:max(lengths)]
    }
    
    # single gaussian (lin)
    funclin1 = function(x, m, s, n){
        f = (n / (s*sqrt(2*pi))) * exp(-((x-m)^2)/(2*(s^2)))
        return(f)
    }
    
    # double gaussian (lin)
    funclin2 = function(x, m1, s1, n1, m2, s2, n2){
        f = ((n1 / (s1*sqrt(2*pi))) * exp(-((x-m1)^2)/(2*(s1^2)))) + ((n2 / (s2*sqrt(2*pi))) * exp(-((x-m2)^2)/(2*(s2^2))))
        return(f)
    }
    
    # single gaussian (log)
    funclog1 = function(x, m, s, n){
        f = n * ((log10(exp(1)))/(x*s*sqrt(2*pi))) * (exp(-((log10(x)-m)^2)/(2*(s^2))))
        return(f)
    }
    
    # double gaussian (log)
    funclog2 = function(x, m1, s1, n1, m2, s2, n2){
        f = (n1 * ((log10(exp(1)))/(x*s1*sqrt(2*pi))) * (exp(-((log10(x)-m1)^2)/(2*(s1^2))))) + (n2 * ((log10(exp(1)))/(x*s2*sqrt(2*pi))) * (exp(-((log10(x)-m2)^2)/(2*(s2^2)))))
        return(f)
    }
    
    # single gaussian (ln)
    funcln1 = function(x, m, s, n){
        f = n * (1/(x*s*sqrt(2*pi))) * (exp(-((log(x)-m)^2)/(2*(s^2))))
        return(f)
    }
    
    # double gaussian (ln)
    funcln2 = function(x, m1, s1, n1, m2, s2, n2){
        f = (n1 * (1/(x*s1*sqrt(2*pi))) * (exp(-((log(x)-m1)^2)/(2*(s1^2))))) + (n2 * (1/(x*s2*sqrt(2*pi))) * (exp(-((log(x)-m2)^2)/(2*(s2^2)))))
        return(f)
    }
    
    # integration steps
    stepslo = x-bw/2
    stepshi = x+bw/2
    
    # integrate
    phi = {}
    for(i in 1:length(x)){
        
        if(ftype=="log"){
            # logged
            if(max(lengths)==1){
                temp = integrate(funclog1, lower=stepslo[i], upper=stepshi[i], m=mean[1], s=sigma[1], n=norm[1], ...)
            }else{
                temp = integrate(funclog2, lower=stepslo[i], upper=stepshi[i], m1=mean[1], s1=sigma[1], n1=norm[1], m2=mean[2], s2=sigma[2], n2=norm[2], ...)
            }
        }else if(ftype=="ln"){
            # natural log
            if(max(lengths)==1){
                temp = integrate(funcln1, lower=stepslo[i], upper=stepshi[i], m=mean[1], s=sigma[1], n=norm[1], ...)
            }else{
                temp = integrate(funcln2, lower=stepslo[i], upper=stepshi[i], m1=mean[1], s1=sigma[1], n1=norm[1], m2=mean[2], s2=sigma[2], n2=norm[2], ...)
            }
        }else{
            # linear
            if(max(lengths)==1){
                temp = integrate(funclin1, lower=stepslo[i], upper=stepshi[i], m=mean[1], s=sigma[1], n=norm[1], ...)
            }else{
                temp = integrate(funclin2, lower=stepslo[i], upper=stepshi[i], m1=mean[1], s1=sigma[1], n1=norm[1], m2=mean[2], s2=sigma[2], n2=norm[2], ...)
            }
        }
        phi = c(phi, temp$value)
        
    }
    
    # return data
    return(phi/bw)
    
}

