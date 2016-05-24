# n=3000*2; followup.length=2; incidence.density=0.01; age.sim="continuous"; random.censoring.rate=0.05; seed=1
sim.dat.tvarying.two = function (
    n, followup.length, incidence.density, 
    age.sim=c("tvaryinggroup","baselinegroup","continuous","bt"),
    random.censoring.rate=0.05, 
    seed) {
    
    set.seed(seed)    
    age.sim <- match.arg(age.sim)    
    
    # initial data frame
    dat=data.frame(
        ptid=1:n,  # a total of n subjects
        trt=rep(c(0,1), each=n/2),  # 1:1 randomization ratio    
        for.non.tvarying.ana=TRUE) # this is added to make it easier to create a new dataset for non-time varying analysis    
    
    # uniform censoring distribution
    if(random.censoring.rate==0) dat$C=followup.length else dat$C=pmin(followup.length, runif(n, 0, followup.length/random.censoring.rate))
    
    # baseline age is almost uniformly distributed between 2.0 and 15.0 in CYD14, in this simulation I change the upper bound to 16 to have more subjects in the oldest age group
    dat$baseline.age=round(runif(n,2,12),3) # round to 3rd digit b/c 1yr=365days, it is necessary to ensure later on that tstart<tstop
    # [0,6)   [6,100)
    dat$agegrp=cut(dat$baseline.age, breaks=c(0,6,100), right=FALSE)
    table(dat$agegrp)
    dat$baseline.agegrp=dat$agegrp
    
    # simulate outcome using baseline variables
    if (age.sim=="continuous") {
        # colnames(design.mat): #  "trt"    "baseline.age"     "trt:baseline.age"   
        design.mat=model.matrix(~trt*baseline.age, dat)[,-1] 
        coef.=c(-0.3, -0.1, -0.08) # determined from simulation results of categorical age group models
    } else {
        # colnames(design.mat): #  "trt"    "agegrp[6,100)"  "trt:agegrp[6,100)"
        design.mat=model.matrix(~trt*agegrp, dat)[,-1] # cannot add -1+ to formula because of the factor
        coef.=c(-0.5, -1, -0.5) # average of cyd14 and koltata
    }    
    exp.linear.predictors = exp(drop(design.mat %*% coef.)) # need to drop b/c linear.predictors is often obtained from matrix product
    baseline.hazard = incidence.density/(0.31-0.02*(followup.length-2)) # this is determined empirically
    ft = -log(1-runif(n)) / (exp.linear.predictors*baseline.hazard) # transform a unif rv to the desired rv
    dat$tstart=0 #we first create a dataset for time varying analysis
    dat$tstop=pmin(ft, dat$C)
    dat$d=ifelse(dat$C>ft, 1, 0) # d = 1 means an event
    
    #print(mean(dat$d))
    
    # create datasets for time varying analysis
    
    # if a 2~5 year old has a followup time that is 6+, we need to make two copies of the subject
    tmp.1=dat[dat$agegrp=="[0,6)" & dat$tstop+dat$baseline.age>=6,]
    nrow(tmp.1)
    if (nrow(tmp.1)>0) {
        # create a copy first
        tmp.2=tmp.1
        # The first copy should always have d=0 
        # the first copy is not used in non-tvaring analysis
        tmp.1$for.non.tvarying.ana=FALSE
        tmp.1$d=0
        tmp.1$tstop=5.9999-tmp.1$baseline.age # 5.9999 instead of 6 because for coxph, each row is an interval of observation (start, stop]
        # the second copy has new agegrp
        tmp.2$tstart=tmp.1$tstop
        tmp.2$agegrp="[6,100)"
        tmp.2$agegrp=factor(tmp.2$agegrp, levels=c("[0,6)","[6,100)"))    
        if (age.sim=="tvaryinggroup" | age.sim=="bt") {
            # the followup time for the second copy need to be generated again using the new predictor
            if(age.sim=="tvaryinggroup"){
                design.mat=model.matrix(~trt*agegrp, tmp.2)[,-1] # cannot add -1 to formula because of the factor
            } else if (age.sim=="bt") {
                design.mat=model.matrix(~trt+agegrp+trt:baseline.agegrp, tmp.2)[,-1] # cannot add -1 to formula because of the factor
            }
            exp.linear.predictors = exp(drop(design.mat %*% coef.)) # need "drop" b/c linear.predictors is often obtained from matrix product
            ft = -log(1-runif(nrow(tmp.2))) / (exp.linear.predictors*baseline.hazard) # transform a unif rv to the desired rv
            tmp.2$tstop = pmin(tmp.2$tstart+ft, tmp.2$C) # no need to re-generate C
            tmp.2$d=ifelse(tmp.2$C>tmp.2$tstop, 1, 0) # d = 1 means an event
        }
        # replace the old copy with the updated copy and the second copy
        tmp.3=rbind(tmp.1, tmp.2)
        tmp.3=tmp.3[order(tmp.3$ptid),]
        dat=rbind(subset(dat, !(dat$agegrp=="[0,6)" & dat$tstop+dat$baseline.age>=6)), tmp.3)
    }
        
    dat$X=dat$tstop    
    invisible (dat)
}
