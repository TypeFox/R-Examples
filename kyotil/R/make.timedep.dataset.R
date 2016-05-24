# X="X"; d="d"; baseline.ageyrs="baseline.age"; t.1=6
make.timedep.dataset = function(dat, X, d, baseline.ageyrs, t.1, t.2=NULL) {
    
    stopifnot(t.1<t.2)
    
    # round to 3rd digit b/c 1yr=365days, it is necessary to ensure later on that tstart<tstop
    dat[[baseline.ageyrs]]=round(dat[[baseline.ageyrs]],3)
    dat$tstart = 0
    dat$tstop = dat[[X]]
    
    # cut into two or three groups depending on t.2
    if (is.null(t.2)) breaks=c(0,t.1,100) else breaks=c(0,t.1,t.2,100)
    dat$.timedep.agegrp=cut(dat[[baseline.ageyrs]], breaks=breaks, right=FALSE)
    dat$.baseline.agegrp=dat$.timedep.agegrp
    .levels=levels(dat$.timedep.agegrp)        
        
    # subset.2 contains subjects entering the study in the young group, but finishes the study in an older group
    # create two copies for each row in subset.2
    age.at.X=dat[[baseline.ageyrs]] + dat[[X]] # this way to calculate age is an approximation, the more accurate calculation can be used, but it is quite time consuming
    subset.2 = dat[[baseline.ageyrs]] < t.1   &   age.at.X>=t.1
    dat.young = subset(dat, subset.2)
    dat.young.cpy=dat.young
    if (nrow(dat.young)>0) {
        # The first copy should have d=0 even for cases
        dat.young[[d]]=0
        dat.young$tstop=(t.1 - 1e-4) - dat.young[[baseline.ageyrs]]  # -1e-4 is necessary b/c in coxph, the format is (start, stop]
        # The second copy should have a new age group
        dat.young.cpy$tstart=dat.young$tstop
        dat.young.cpy$.timedep.agegrp=.levels[2]
        dat.young.cpy$.timedep.agegrp=factor(dat.young.cpy$.timedep.agegrp, levels=.levels)    
    }
    
    if (is.null(t.2)) {
        out=rbind(
            subset(dat, !subset.2), 
            dat.young, dat.young.cpy
        )
    
    } else {
        # subset.1 contains subjects entering the study in the middle group, but finishes the study in the old group
        # create two copies for each row in subset.1
        subset.1 = dat[[baseline.ageyrs]]<t.2 & dat[[baseline.ageyrs]]>=t.1   &   age.at.X>=t.2
        dat.middle = subset(dat, subset.1)
        dat.middle.cpy=dat.middle
        if(nrow(dat.middle)>0) {
            # The first copy should have d=0 even for cases
            dat.middle[[d]]=0
            dat.middle$tstop=(t.2 - 1e-4) - dat.middle[[baseline.ageyrs]]  # -1e-4 is necessary b/c in coxph, the format is (start, stop]
            # The second copy should have a new age group
            dat.middle.cpy$tstart=dat.middle$tstop
            dat.middle.cpy$.timedep.agegrp=.levels[3]
            dat.middle.cpy$.timedep.agegrp=factor(dat.middle.cpy$.timedep.agegrp, levels=.levels)    
        }
    
        out=rbind(
            subset(dat, !(subset.1 | subset.2)), 
            dat.middle, dat.middle.cpy,
            dat.young, dat.young.cpy
        )
    }
    
}
