brkdn.plot<-function(vars,groups=NA,obs=NA,data,mct="mean",md="std.error",
 stagger=NA,dispbar=TRUE,main="Breakdown plot",xlab=NA,ylab=NA,xaxlab=NA,
 ylim=NA,type="b",pch=1,lty=1,col=par("fg"),staxx=FALSE,...) {

 if(class(vars) == "formula") {
  formbits<-all.vars(vars)
  vars<-formbits[1]
  groups<-formbits[2]
  obs<-formbits[3]
 }
 if(is.na(obs)) {
  if(is.na(groups[1]))
   stop("Must have at least one factor to subset data")
  bygroup<-as.factor(data[[groups]])
  grouplevels<-levels(bygroup)
  ngroups<-length(grouplevels)
  nobs<-length(vars)
  obs.pos<-1:nobs
  obslevels<-1:nobs
 }
 else {
  if(is.numeric(data[[obs]])) {
   obs.pos<-obslevels<-sort(unique(data[[obs]]))
   nobs<-length(obslevels)
  }
  else {
   byobs<-as.factor(data[[obs]])
   obslevels<-levels(byobs)
   nobs<-length(obslevels)
   obs.pos<-1:nobs
  }
  if(is.na(groups)) {
   ngroups<-length(vars)
   grouplevels<-1:ngroups
  }
  else {
   bygroup<-as.factor(data[[groups]])
   grouplevels<-levels(bygroup)
   ngroups<-length(grouplevels)
   if(length(vars) > 1) {
    warning("Group and observation factors are present, only vars[1] is plotted")
    vars<-vars[1]
   }
  }
 }
 brkdn<-list(matrix(NA,nrow=ngroups,ncol=nobs),
  matrix(NA,nrow=ngroups,ncol=nobs))
 if(is.na(groups)) {
  if(is.na(xlab)) xlab<-"Observation"
  xat<-1:nobs
  if(is.na(xaxlab[1])) xaxlab<-obslevels
  for(group in 1:ngroups) {
   for(ob in 1:nobs) {
    thisbit<-data[[vars[group]]][data[[obs]] == obslevels[ob]]
    if(length(thisbit)) {
     if(length(thisbit) > 1) {
      brkdn[[1]][group,ob]<-do.call(mct,list(thisbit,na.rm=TRUE))
      if(!is.na(md))
       brkdn[[2]][group,ob]<-do.call(md,list(thisbit,na.rm=TRUE))
     }
     else brkdn[[1]][group,ob]<-thisbit
    }
   }
  }
 }
 else {
  if(is.na(obs)) {
   if(is.na(xlab)) xlab<-"Variable"
   xat<-1:length(vars)
   if(is.na(xaxlab[1])) xaxlab<-vars
   for(group in 1:ngroups) {
    for(ob in 1:nobs) {
     thisbit<-data[[vars[ob]]][data[[groups]] == grouplevels[group]]
     if(length(thisbit)) {
      if(length(thisbit) > 1) {
       brkdn[[1]][group,ob]<-do.call(mct,list(thisbit,na.rm=TRUE))
       if(!is.na(md))
       brkdn[[2]][group,ob]<-do.call(md,list(thisbit,na.rm=TRUE))
      }
      else brkdn[[1]][group,ob]<-thisbit
     }
    }
   }
  }
  else {
   if(is.na(xlab)) xlab<-"Observation"
   xat<-obs.pos
   if(is.na(xaxlab[1])) xaxlab<-obslevels
   for(group in 1:ngroups) {
    for(ob in 1:nobs) {
     thisbit<-data[[vars]][data[[groups]] == grouplevels[group] &
       data[[obs]] == obslevels[ob]]
     if(length(thisbit)) {
      if(length(thisbit) > 1) {
       brkdn[[1]][group,ob]<-do.call(mct,list(thisbit,na.rm=TRUE))
       if(!is.na(md))
        brkdn[[2]][group,ob]<-do.call(md,list(thisbit,na.rm=TRUE))
      }
      else brkdn[[1]][group,ob]<-thisbit
     }
    }
   }
  }
 }
 if(is.na(ylim[1])) {
  ylim<-range(brkdn[[1]],na.rm=TRUE)
  if(!is.na(md)) {
   dlim<-c(min(brkdn[[1]]-brkdn[[2]],na.rm=TRUE),
     max(brkdn[[1]]+brkdn[[2]],na.rm=TRUE))
   ylim<-c(min(c(ylim[1],dlim[1])),max(c(ylim[2],dlim[2])))
  }
 }
 groupdiv<-ifelse(ngroups < 3,1,ngroups-2)
 if(is.na(stagger)) stagger<-0.025-groupdiv*0.0025
 if(is.na(ylab)) {
  if(length(vars) == 1) ylab<-vars[1]
  else ylab<-paste(vars,collapse=" and ")
 }
 plot(0,xlim=c(obs.pos[1]-0.5,obs.pos[nobs]+0.5),main=main,
  xlab=xlab,ylab=ylab,ylim=ylim,type="n",axes=FALSE,...)
 box()
 if(staxx) staxlab(at=xat,labels=xaxlab)
 else axis(1,at=xat,labels=xaxlab)
 axis(2)
 if(length(pch) < ngroups) pch<-rep(pch,length.out=ngroups)
 if(length(col) < ngroups) col<-rep(col,length.out=ngroups)
 if(length(lty) < ngroups) lty<-rep(lty,length.out=ngroups)
 offinc<-stagger*diff(par("usr")[c(1,2)])
 offset<-0
 arrow.cap<-0.01-(groupdiv*0.001)
 for(group in 1:ngroups) {
  bg<-col[group]
  if(pch[group] > 20 && pch[group] < 26) col[group]<-par("fg")
  points(obs.pos+offset,brkdn[[1]][group,],type=type,col=col[group],
   pch=pch[group],lty=lty[group],bg=bg,...)
  if(dispbar)
   dispersion(obs.pos+offset,brkdn[[1]][group,],brkdn[[2]][group,],
    arrow.cap=arrow.cap,col=col[group],...)
  offset<-ifelse(offset<0,-offset,-offset-offinc)
 }
 names(brkdn)<-c(mct,md)
 return(brkdn)
}
