combrad = function(mag, re, n){
    # merge multiple sersic functions to calculate a global flux radius
    # base function
    calccombrad = function(mag, re, n, radii){
        maxlength = max(length(mag), length(re), length(n))[1]
        proflum = rep(0, length(radii))
        proflumtot = 0
        profs = {}
        for(i in 1:maxlength){
            prof = sersic(mag=mag[i], re=re[i], n=n[i], e=0, r=radii)
            profs = c(profs, list(prof))
            proflum = proflum + prof$lum
            proflumtot = proflumtot + prof$lumtot
        }
        frac = proflum/proflumtot
        hlrad = radii[which.min(abs(frac-0.5))]
        return(hlrad)
    }
    # function to loop over base function, for accuracy
    loopcombrad = function(mag, re, n, zooms=c(0.1,0.05,0.01)){
        refradius = max(re)[1]
        radii = refradius * seq(0,20,len=(1E3)+1)
        hlbest = calccombrad(mag, re, n, radii)
        if(length(hlbest)!=0){
            for(i in 1:length(zooms)){
                lo = hlbest * (1-zooms[i])
                hi = hlbest * (1+zooms[i])
                radii = seq(lo,hi,len=(1E3)+1)
                hlbest = calccombrad(mag, re, n, radii)
            }
        }else{
            hlbest = NA
        }
        return(hlbest)
    }
    # fix data lists of differing lengths
    data = list(mag, re, n)
    ncomp = c(length(mag), length(re), length(n))
    ncompmax = max(ncomp)
    if(any(ncomp != ncompmax)){
        bad = which(ncomp != ncompmax)
        for(i in 1:length(bad)){
            badlen = ncomp[bad[i]]
            goodbaddiff = ncompmax - badlen
            added = {}
            for(j in 1:goodbaddiff){
                added = c(added, list(NA))
            }
            data[[bad[[i]]]] = c(data[[bad[[i]]]], added)
        }
    }
    # calculate input data lengths
    datalen = {}
    for(i in 1:3){
        for(j in 1:ncompmax){
            datalen = c(datalen, length(data[[i]][[j]]))
        }
    }
    datalenmax = max(datalen)
    # loop over each component
    inmag = {}
    inre = {}
    inindex = {}
    for(i in 1:ncompmax){
        mm = data[[1]][[i]]
        rr = data[[2]][[i]]
        nn = data[[3]][[i]]
        if(length(mm) != datalenmax){mm = c(mm,rep(NA,(datalenmax-length(mm))))}
        if(length(rr) != datalenmax){rr = c(rr,rep(NA,(datalenmax-length(rr))))}
        if(length(nn) != datalenmax){nn = c(nn,rep(NA,(datalenmax-length(nn))))}
        inmag = rbind(inmag, mm)
        inre = rbind(inre, rr)
        inindex = rbind(inindex, nn)
    }
    # calculate half-light radii for each object
    hlrads = {}
    for(i in 1:datalenmax){
        h = loopcombrad(mag=inmag[,i], re=inre[,i], n=inindex[,i])
        hlrads = c(hlrads, h)
    }
    return(hlrads)
}

