# calculate number (total), number density (phi) and poissonian errors
schechter.bin = function(data, vmax = NA, range = range(data), lim1 = NA, lim2 = NA, numlim = 1, volume = max(vmax), bw = 0.1, null = 1E-9){
    
    # data vectors
    bins = seq(range[1], range[2], by=bw)
    binmid = bins[-1]-(bw/2)
    binlo = bins[1:(length(bins)-1)]
    binhi = bins[2:length(bins)]
    
    # volume
    if(is.na(volume[1])){volume = 1}
    if(!is.na(vmax[1])){
        vmin = unlist(lapply(vmax, min, volume))
        vweight = volume / vmin
    }else{
        vweight = rep(1, len=length(data))
    }
    
    # calculate bin sample statistics
    num = {}
    den = {}
    err = {}
    for(i in 1:length(binmid)){
        
        # sub-sample
        samp = vweight[data>binlo[i] & data<binhi[i]]
        
        # calculate bin values
        sampnum = sum(samp)
        sampden = (sampnum/volume)/bw
        samperr = (sqrt(sampnum)/volume)/bw
        
        # add to data vectors
        num = c(num,sampnum)
        den = c(den,sampden)
        err = c(err,samperr)
        
    }
    
    # calculate upper and lower errors
    errlo = den-err
    errhi = den+err
    if(any(errlo==0)){
        errlo[which(errlo==0)] = null
    }
    if(any(errhi==0)){
        errhi[which(errhi==0)] = null
    }
    
    # schechter (to be) fit values
    if(any(num<=numlim)){
        bad = which(num<=numlim)
        fitbinmid = binmid[-bad]
        fitbinlo = binlo[-bad]
        fitbinhi = binhi[-bad]
        fitnum = num[-bad]
        fitden = den[-bad]
        fiterr = err[-bad]
        fiterrlo = errlo[-bad]
        fiterrhi = errhi[-bad]
    }else{
        fitbinmid = binmid
        fitbinlo = binlo
        fitbinhi = binhi
        fitnum = num
        fitden = den
        fiterr = err
        fiterrlo = errlo
        fiterrhi = errhi
    }
    
    # impose any upper and lower limits
    bad = {}
    if(!is.na(lim1)){
        if(any(fitbinlo<lim1)){
            bad = c(bad,which(fitbinlo<lim1))
        }
    }
    if(!is.na(lim2)){
        if(any(fitbinhi>lim2)){
            bad = c(bad,which(fitbinhi>lim2))
        }
    }
    if(length(bad)>0){
        fitbinmid = fitbinmid[-bad]
        fitbinlo = fitbinlo[-bad]
        fitbinhi = fitbinhi[-bad]
        fitnum = fitnum[-bad]
        fitden = fitden[-bad]
        fiterr = fiterr[-bad]
        fiterrlo = fiterrlo[-bad]
        fiterrhi = fiterrhi[-bad]
    }
    
    # return results
    return(list(bins=bins, binmid=binmid, binlo=binlo, binhi=binhi, num=num, den=den, err=err, errlo=errlo, errhi=errhi, fitbinmid=fitbinmid, fitbinlo=fitbinlo, fitbinhi=fitbinhi, fitnum=fitnum, fitden=fitden, fiterr=fiterr, fiterrlo=fiterrlo, fiterrhi=fiterrhi))
    
}
