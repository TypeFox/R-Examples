schechter.ellipse = function(data, vmax = NA, knee, slope, norm, chi2, datarange = NA, kneerange = c(-24,-16), sloperange = c(-2,1.5), kneeofflims = NA, slopeofflims = NA, kneestep = 0.5, slopestep = 0.1, kneesteps = NA, slopesteps = NA, lim1 = NA, lim2 = NA, numlim = 1, method = "nlminb", volume = max(vmax), bw = 0.1, mag = FALSE, log = FALSE, null = 1E-9, ...){
    
    #load("../vollim.img"); data = dat[dat[,"R_SERSIC"]>-50,"R_SERSIC"]; knee = c(-20.810520443); slope = c(-1.098815208, 0.657166733); norm = c(0.005030941, 0.002933632); chi2=34.716425374; kneerange = c(-24,-15); sloperange = c(-2,1.5); kneestep=1; slopestep=0.5; lim1=NA; lim2=-17.4; numlim=3; method="nlminb"; volume=224555.3; bw=0.25; mag=TRUE; log=FALSE; null = 1E-9
    
    #load("../vollim.img"); data = dat[dat[,"R_SERSIC"]>-50,"R_SERSIC"]; knee = -21.61014746; slope = -1.10482398; norm = 0.00436335; chi2=43.67826350; datarange=c(-24,-15); kneerange = c(-24,-16); sloperange = c(-2,1.5); kneeofflims = c(-3,1); slopeofflims = c(-1,1.1); kneestep=1; slopestep=0.5; kneesteps = 20; slopesteps = 20; lim1=NA; lim2=-17.4; numlim=3; method="nlminb"; volume=224555.3; bw=0.25; mag=TRUE; log=FALSE; null = 1E-9
    
#    # check inputs have same length (for double schechter fits)
#    if(any(is.na(knee))){knee = knee[-which(is.na(knee))]}
#    if(any(is.na(slope))){slope = slope[-which(is.na(slope))]}
#    if(any(is.na(norm))){norm = norm[-which(is.na(norm))]}
#    lengths = c(length(knee), length(slope), length(norm))
#    if(max(lengths)!=min(lengths)){
#        knee = rep(knee,max(lengths))[1:max(lengths)]
#        slope = rep(slope,max(lengths))[1:max(lengths)]
#        norm = rep(norm,max(lengths))[1:max(lengths)]
#    }
    
    # volume
    if(is.na(volume[1])){volume = 1}
    
    # data range check
    if(is.na(datarange[1])){
        datarange = kneerange
    }
    
    # mag check
    if(mag){
        log = FALSE
    }
    
    # set up knee/slope ranges
    if(!is.na(kneeofflims[1])){
        kneerange = c( (min(knee)+min(kneeofflims)), (max(knee)+max(kneeofflims)) )
    }
    if(!is.na(slopeofflims[1])){
        sloperange = c( (min(slope)+min(slopeofflims)), (max(slope)+max(slopeofflims)) )
    }
    
    # set up knee/slope intervals
    if(!is.na(kneesteps[1])){
        knees = seq(sort(kneerange)[1], sort(kneerange)[2], len=kneesteps[1])
    }else{
        knees = seq(sort(kneerange)[1], sort(kneerange)[2], by=abs(kneestep))
    }
    if(!is.na(slopesteps[1])){
        slopes = seq(sort(sloperange)[1], sort(sloperange)[2], len=slopesteps[1])
    }else{
        slopes = seq(sort(sloperange)[1], sort(sloperange)[2], by=abs(slopestep))
    }
    
    # set up fitting grid
    fullgrid = expand.grid(knees,slopes)
    
    # grid
    grid = fullgrid
    
    # run error fitting subroutine
    res = .schechter.ellipse.fit(grid=grid, vmax=vmax, data=data, knee=knee, slope=slope, norm=norm, datarange=datarange, lim1=lim1, lim2=lim2, numlim=numlim, method=method, volume=volume, bw=bw, mag=mag, log=log, null=null)
    
    # construct res1 and res2 (if applicable) error matrices
    res1 = t(matrix(res$res1, ncol=length(knees), nrow=length(slopes), byrow=TRUE))
    if(is.null(res$res2)){
        res2 = res$res2 
    }else{
        res2 = t(matrix(res$res2, ncol=length(knees), nrow=length(slopes), byrow=TRUE))
    }
    
    # define sigma limits - See numerical recipies and also Lin et al, 1996 for a good summary
    # note only 2 DoF used here, phi* must be determined in some other fashion
    s1 = as.numeric(chi2)+2.30
    s2 = as.numeric(chi2)+6.17
    s3 = as.numeric(chi2)+11.8
    
    #contour(knees, slopes, res1, levels=c(s1,s2,s3)); abline(v=knee[1]); abline(h=slope[1])
    
    # calculate error bounds
    s1r1cols = which((res1-s1)<=0, arr.ind=TRUE)
    s2r1cols = which((res1-s2)<=0, arr.ind=TRUE)
    s3r1cols = which((res1-s3)<=0, arr.ind=TRUE)
    if(length(s1r1cols[,1])==0){s1r1cols = rbind(s1r1cols,c(NA,NA))}
    if(length(s2r1cols[,1])==0){s2r1cols = rbind(s2r1cols,c(NA,NA))}
    if(length(s3r1cols[,1])==0){s3r1cols = rbind(s3r1cols,c(NA,NA))}
    kneelo = suppressWarnings(c(min(s1r1cols[,1]), min(s2r1cols[,1]), min(s3r1cols[,1])))
    kneehi = suppressWarnings(c(max(s1r1cols[,1]), max(s2r1cols[,1]), max(s3r1cols[,1])))
    slopelo = suppressWarnings(c(min(s1r1cols[,2]), min(s2r1cols[,2]), min(s3r1cols[,2])))
    slopehi = suppressWarnings(c(max(s1r1cols[,2]), max(s2r1cols[,2]), max(s3r1cols[,2])))
    kneelo1 = suppressWarnings(knees[kneelo]-knee[1])
    kneehi1 = suppressWarnings(knees[kneehi]-knee[1])
    slopelo1 = suppressWarnings(slopes[slopelo]-slope[1])
    slopehi1 = suppressWarnings(slopes[slopehi]-slope[1])
    
    if(!is.null(res2[1])){
        s1r2cols = which((res2-s1)<=0, arr.ind=TRUE)
        s2r2cols = which((res2-s2)<=0, arr.ind=TRUE)
        s3r2cols = which((res2-s3)<=0, arr.ind=TRUE)
        if(length(s1r2cols[,1])==0){s1r2cols = rbind(s1r2cols,c(NA,NA))}
        if(length(s2r2cols[,1])==0){s2r2cols = rbind(s2r2cols,c(NA,NA))}
        if(length(s3r2cols[,1])==0){s3r2cols = rbind(s3r2cols,c(NA,NA))}
        kneelo = suppressWarnings(c(min(s1r2cols[,1]), min(s2r2cols[,1]), min(s3r2cols[,1])))
        kneehi = suppressWarnings(c(max(s1r2cols[,1]), max(s2r2cols[,1]), max(s3r2cols[,1])))
        slopelo = suppressWarnings(c(min(s1r2cols[,2]), min(s2r2cols[,2]), min(s3r2cols[,2])))
        slopehi = suppressWarnings(c(max(s1r2cols[,2]), max(s2r2cols[,2]), max(s3r2cols[,2])))
        kneelo2 = suppressWarnings(knees[kneelo]-knee[length(knee)])
        kneehi2 = suppressWarnings(knees[kneehi]-knee[length(knee)])
        slopelo2 = suppressWarnings(slopes[slopelo]-slope[length(slope)])
        slopehi2 = suppressWarnings(slopes[slopehi]-slope[length(slope)])
    }else{
        kneelo2 = kneehi2 = slopelo2 = slopehi2 = NA
    }
        
    # return error data
    return(list(knees=knees, slopes=slopes, res1=res1, res2=res2, s1=s1, s2=s2, s3=s3, kneelo1=kneelo1, kneehi1=kneehi1, slopelo1=slopelo1, slopehi1=slopehi1, kneelo2=kneelo2, kneehi2=kneehi2, slopelo2=slopelo2, slopehi2=slopehi2))
    
}

.schechter.ellipse.fit = function(grid, vmax, data, knee, slope, norm, datarange, lim1, lim2, numlim, method, volume, bw, mag, log, null){
    
    # inputs lengths
    lengths = c(length(knee), length(slope), length(norm))
    
    i = 1
    
    # loop over each grid points
    res1 = {}
    res2 = {}
    for(i in 1:length(grid[,1])){
        
        # setup
        cat("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",i,"/",length(grid[,1]),"")
        vars = grid[i,]
        
        j = 1
        
        # loop over each input schechter function
        for(j in 1:max(lengths)){
            
            # if first schechter function (or single)
            if(j==1){
                
                # setup inputs
                k = knee
                s = slope
                k[1] = vars[,1]
                s[1] = vars[,2]
                
                # fit function
                fit = schechter.fit(data=data, vmax=vmax, knee=k, slope=s, norm=norm, range=datarange, lim1=lim1, lim2=lim2, numlim=numlim, method=method, volume=volume, bw=bw, mag=mag, log=log, null=null, fixk1=TRUE, fixs1=TRUE, error="none")
                
                # add chi2 to results
                res1 = c(res1, fit$chi2)
                
            }else{
                
                # setup inputs
                k = knee
                s = slope
                k[length(knee)] = vars[,1]
                s[length(knee)] = vars[,2]
                
                # fit function
                if(length(knee)==1){
                    fit = schechter.fit(data=data, vmax=vmax, knee=k, slope=s, norm=norm, range=datarange, lim1=lim1, lim2=lim2, numlim=numlim, method=method, volume=volume, bw=bw, mag=mag, log=log, null=null, fixk1=TRUE, fixs2=TRUE, error="none")
                }else{
                    fit = schechter.fit(data=data, vmax=vmax, knee=k, slope=s, norm=norm, range=datarange, lim1=lim1, lim2=lim2, numlim=numlim, method=method, volume=volume, bw=bw, mag=mag, log=log, null=null, fixk2=TRUE, fixs2=TRUE, error="none")
                }
                
                # add chi2 to results
                res2 = c(res2, fit$chi2)
                
            }
            
        }
        
    }
    
    # return error data
    return(list(grid=grid, res1=res1, res2=res2))
    
}

