phase.sync <- function (t1, t2, nrands = 0, mod = 1, nbreaks = 10, 
                        mins = FALSE, quiet = FALSE) {
  
  if (NCOL(t1)==1 | NCOL(t2)==1) {
    t1=cbind(1:NROW(t1), t1)
    t2=cbind(1:NROW(t2), t2)    
  }
  p=phase.sync.aux(t1, t2, mins=mins)
  
  if (nrands == 0) {
    results=p
  }
  else {
    rands=numeric(length=nrands+1)*NA
    if (mod == 1) {
      column=3
      breaks=seq(from=0, to=2*pi, length.out=10)
    }
    else {
      column=4
      breaks=seq(from=-pi, to=pi, length.out=10)
    }
    h.obs=hist(p$deltaphase[, column], breaks=breaks, plot=FALSE)
    p.obs=h.obs$counts/sum(h.obs$counts)
    nbins.obs=length(h.obs$counts)
    S.obs=-sum(p.obs*log(p.obs), na.rm=TRUE)
    Smax.obs=log(nbins.obs)
    Q.obs=(Smax.obs-S.obs)/Smax.obs
    
    # Determine transition probabilities
    distr.t1 <- cut(t1[,2], quantile(t1[,2], seq(0, 1, len = (nbreaks+1))), 
                    include.lowest = TRUE, labels=FALSE)
    distr.t2 <- cut(t2[,2], quantile(t2[,2], seq(0, 1, len = (nbreaks+1))), 
                    include.lowest = TRUE, labels=FALSE)  
    trans.t1=matrix(nrow=nbreaks, ncol=nbreaks, 0)
    trans.t2=matrix(nrow=nbreaks, ncol=nbreaks, 0)
    
    for (i in 1:(NROW(t1)-1)) {
      trans.t1[distr.t1[i], distr.t1[i+1]]=trans.t1[distr.t1[i], distr.t1[i+1]]+1      
      trans.t2[distr.t2[i], distr.t2[i+1]]=trans.t2[distr.t2[i], distr.t2[i+1]]+1      
    }
    trans.t1=trans.t1/rowSums(trans.t1)  
    trans.t2=trans.t2/rowSums(trans.t2)  
    if (!quiet)
      prog.bar=txtProgressBar(min = 0, max = nrands, style = 3)
    for (r in 1:nrands) {
      # Surrogate randomization (Cazelles and Stone 2003)
      surr.t1=surrogate.ts(ts=t1, distr.ts=distr.t1, nbreaks=nbreaks)$surr.ts
      surr.t2=surrogate.ts(ts=t2, distr.ts=distr.t2, nbreaks=nbreaks)$surr.ts
      p.rand=phase.sync(surr.t1, surr.t2)
      rand.h=hist(p.rand$deltaphase[, column], breaks=breaks, plot=FALSE)
      rand.p=rand.h$counts/sum(rand.h$counts, na.rm=TRUE)
      rand.nbins=length(rand.h$counts)  
      rand.S=-sum(rand.p*log(rand.p), na.rm=TRUE)
      rand.Smax=log(rand.nbins)
      rands[r]=(rand.Smax-rand.S)/rand.Smax
      if (!quiet)
        setTxtProgressBar(prog.bar, r)
    }
    rands[r+1]=Q.obs
    pValue = sum (rands >= Q.obs)/(nrands+1)  
    # ICDF
    o=sort(rands)
    icdf=data.frame(Q=o, icdf=sapply(o, FUN=function (x) {sum(rands >= x)/(nrands+1)}))
    results=list(Q.obs=Q.obs, pval=pValue, rands=rands, phases1=p$phases1,
                 phases2=p$phases2, deltaphase=p$deltaphase, icdf=icdf)
  }
  class(results)="phase"
  return (results)
}

phase.sync.aux <- function (t1, t2, mins = FALSE) {
  # Find the min/max in both timeseries
  min.max1=find.minmax(t1)
  min.max2=find.minmax(t2)
  # timesteps, densities, phase
  phases1=matrix(nrow=NROW(t1), ncol=3, NA)
  phases2=matrix(nrow=NROW(t2), ncol=3, NA)
  phases1[,1:2]=t1
  phases2[,1:2]=t2
  if (mins) {
    v1=min.max1$mins
    v2=min.max2$mins
  }
  else {
    v1=min.max1$maxs
    v2=min.max2$maxs
  }
  # Locations of mins/maxs
  locs1=v1$index
  locs2=v2$index
  ## Range of values over which to interpolate
  range1=locs1[1]:locs1[length(locs1)]
  range2=locs2[1]:locs2[length(locs2)]
  # Assign phase values to mins/maxs
  phases1[locs1, 2:3]=cbind(v1[, 2], 
                            seq(from=0, by=2*pi, to=(NROW(v1)-1)*2*pi))
  phases2[locs2, 2:3]=cbind(v2[, 2], 
                            seq(from=0, by=2*pi, to=(NROW(v2)-1)*2*pi))
  # Interpolate phase values between successive mins/maxs
  phases1[range1, 3]=
    approx(x=phases1[, 1], y=phases1[, 3], 
           n=length(range1))$y
  phases2[range2, 3]=
    approx(x=phases2[,1], y=phases2[,3],
           n=length(range2))$y
  
  deltaphase=matrix(nrow=NROW(t1), ncol=4, NA)
  phase_diff=phases1[,3]-phases2[,3]
  mod_phase_diff1=phase_diff %% (2*pi)
  mod_phase_diff2=((phase_diff+pi) %% (2*pi))-pi
  deltaphase=cbind(1:NROW(t1), phase_diff, mod_phase_diff1, mod_phase_diff2)
  colnames(deltaphase)=c("timestep", "phasediff", "mod_phase_diff_2pi", 
                         "mod_phase_diff_pi")
  colnames(phases1)=c("timestep", "val", "phase")
  colnames(phases2)=c("timestep", "val", "phase")
  return (list(phases1=as.data.frame(phases1), 
               phases2=as.data.frame(phases2), 
               deltaphase=as.data.frame(deltaphase)))
}
