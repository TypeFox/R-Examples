lumdens = function(knee, slope, norm, msun = solar("r"), mag = TRUE, log = FALSE){
    
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
    
    # calculate
    if(mag){
        j = norm * (10^(-0.4*(knee-msun))) * gamma(slope+2)
    }else{
        if(log){
            j = norm * (10^knee) * gamma(slope+2)
        }else{
            j = norm * knee * gamma(slope+2)
        }
    }
    
    # return results
    return(sum(j))
    
}

