#Implementations of bootstrap algorithms for the AVE- and MAX-tests (Hung, 2000)
#Metric data: AVE-test based on Student's t-test in bifactorial designs
avestudent2Boot<-function(C,nboot=NULL,simerror=NULL,...){
  if(is.binary(C@data)) warning("For binary data, the binomial test should be used.")
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  cb1<-function(a,b){((a-1)*C@D[2])+b}
  dauer<-proc.time()[3] #Measure duration of simulation
  if(all(C@D==c(1,1))){
    pv<-maxstudent2Boot(C,nboot=nboot,simerror=simerror,...)
    return(new("avetest",pave=pv@pmax,tave=max(pv@tmax),nboot=pv@nboot,simerror=pv@simerror))
  }
  tmin<-numeric(0)
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){
    if(mean(C@data[[cb(a,0)]])>=mean(C@data[[cb(0,b)]])){
      tmin[cb1(a,b)]<-ttest(C@data[[cb(a,b)]],C@data[[cb(a,0)]])
    }
    else{
      tmin[cb1(a,b)]<-ttest(C@data[[cb(a,b)]],C@data[[cb(0,b)]])
    }
  }}
  Y<-carpetdaten(C)
  taveSc<-mean(tmin)/sqrt(var(tmin)/(C@D[1]*C@D[2]))
  results<-.Call("avestudent2",Y,matrix(C@n,nrow=(C@D[1]+1),byrow=TRUE),
                 c(nboot,C@D),c(taveSc,simerror),PACKAGE="bifactorial")
  dauer<-proc.time()[3]-dauer
  new("avetest",stat=mean(tmin),p=round(results[[1]]/results[[2]],4),method="Bootstrap",nboot=results[[2]],
      simerror=sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3)),duration=dauer)
}
#MAX-test based on Student's t-test in bifactorial designs
maxstudent2Boot<-function(C,nboot=NULL,simerror=NULL,...){
  pv<-mintest(C,test="ttest",nboot=nboot,simerror=simerror,...)
  pmax<-round(min(pv@p),4)
  new("maxtest",p=pmax,stat=max(pv@stat),name=pv@gnames[pv@p==min(pv@p)],nboot=pv@nboot,
      simerror=pv@simerror,method="bootstrap",duration=pv@duration)
}
#AVE-test based on Student's t-test in trifactorial designs
avestudent3Boot<-function(C,nboot=NULL,simerror=NULL,...){
  if(is.binary(C@data)) warning("For binary data, the binomial test should be used.")
  cb<-function(a,b,c){(a*(C@D[2]+1)*(C@D[3]+1))+(b*(C@D[3]+1))+c+1}
  cb1<-function(a,b,c){((a-1)*(C@D[2])*(C@D[3]))+((b-1)*(C@D[3]))+c}
  dauer<-proc.time()[3] #Measure duration of simulation
  if(all(C@D==c(1,1,1))){
    pv<-maxstudent3Boot(C,nboot=nboot,simerror=simerror,...)
    return(new("avetest",p=pv@pmax,stat=max(pv@tmax),nboot=pv@nboot,simerror=pv@simerror))
  }
  tmin<-numeric(0)
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){for(c in 1:C@D[3]){
    maxmarg<-max(mean(C@data[[cb(0,b,c)]]),mean(C@data[[cb(a,0,c)]]),mean(C@data[[cb(a,b,0)]]))
    if(maxmarg==mean(C@data[[cb(0,b,c)]])){
      tmin[cb1(a,b,c)]<-ttest(C@data[[cb(a,b,c)]],C@data[[cb(0,b,c)]])
    }
    if(maxmarg==mean(C@data[[cb(a,0,c)]])){
      tmin[cb1(a,b,c)]<-ttest(C@data[[cb(a,b,c)]],C@data[[cb(a,0,c)]])
    }
    if(maxmarg==mean(C@data[[cb(a,b,0)]])){
      tmin[cb1(a,b,c)]<-ttest(C@data[[cb(a,b,c)]],C@data[[cb(a,b,0)]])
    }
  }}}
  Y<-cubedaten(C)
  taveSc<-mean(tmin)/sqrt(var(tmin)/(C@D[1]*C@D[2]*C@D[3]))
  results<-.Call("avestudent3",Y,C@n,c(nboot,C@D),c(taveSc,simerror),PACKAGE="bifactorial")
  dauer<-proc.time()[3]-dauer
  new("avetest",stat=mean(tmin),p=round(results[[1]]/results[[2]],4),
      nboot=results[[2]],simerror=sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3)),
      method="bootstrap",duration=dauer)
}
#MAX-test based on Student's t-test in trifactorial designs
maxstudent3Boot<-function(C,nboot=NULL,simerror=NULL,...){
  pv<-mintest(C,test="ttest",nboot=nboot,simerror=simerror,...)
  pmax<-round(min(pv@p),4)
  new("maxtest",p=pmax,stat=max(pv@stat),name=pv@gnames[pv@p==min(pv@p)],
      nboot=pv@nboot,simerror=pv@simerror,method="bootstrap",duration=pv@duration)
}
#Binary data: AVE-test based on Student's t-test in bifactorial designs
avebinomial2Boot<-function(C,nboot=NULL,simerror=NULL,...){
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  cb1<-function(a,b){((a-1)*C@D[2])+b}
  dauer<-proc.time()[3] #Measure duration of simulation
  if(!is.binary(C@data)) stop("This test needs binary data.")
  n<-matrix(C@n,nrow=C@D[1]+1,byrow=TRUE)
  p<-matrix(nrow=C@D[1]+1,ncol=C@D[2]+1)
  zmin<-matrix(nrow=C@D[1],ncol=C@D[2])
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){p[a+1,b+1]=mean(C@data[[cb(a,b)]])}}
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){
    zmin[a,b]<-if(p[a+1,1]>=p[1,b+1]){ztest(n[a+1,b+1],n[a+1,1],p[a+1,b+1],p[a+1,1])}
    else{ztest(n[a+1,b+1],n[1,b+1],p[a+1,b+1],p[1,b+1])}
  }}
  zave<-mean(zmin)/sqrt(var(as.numeric(zmin))/(C@D[1]*C@D[2]))
  if(is.null(simerror)) simerror<-9 #This is necessary for interfacing to C++
  if(is.null(nboot)) nboot<-900
  results<-.Call("avebinomial2",n,p,c(nboot,C@D),c(zave,simerror),PACKAGE="bifactorial")
  simerror<-sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3))
  dauer<-proc.time()[3]-dauer
  new("avetest",stat=zave,p=round(results[[1]]/results[[2]],4),
      nboot=results[[2]],simerror=simerror,method="Bootstrap",duration=dauer)
}
#MAX-test based on a Z-statistic in bifactorial designs
maxbinomial2Boot<-function(C,nboot=NULL,simerror=NULL,...){
  pv<-mintest(C,test="ztest",nboot=nboot,simerror=simerror,...)
  pmax<-round(min(pv@p),4)
  new("maxtest",p=pmax,stat=max(pv@stat),name=pv@gnames[pv@p==pmax],nboot=pv@nboot,simerror=pv@simerror,method="bootstrap",duration=pv@duration)
}
#AVE-test based on a Z-statistic in trifactorial designs
avebinomial3Boot<-function(C,nboot=NULL,simerror=NULL,...){
  cb<-function(a,b,c){(a*(C@D[2]+1)*(C@D[3]+1))+(b*(C@D[3]+1))+c+1}
  cb1<-function(a,b,c){((a-1)*(C@D[2])*(C@D[3]))+((b-1)*(C@D[3]))+c}
  if(!is.binary(C@data)) stop("This test needs binary data.")
  dauer<-proc.time()[3] #Measure duration of simulation
  p<-zmin<-numeric(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){for(c in 0:C@D[3]){p[cb(a,b,c)]=mean(C@data[[cb(a,b,c)]])}}}
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){for(c in 0:C@D[3]){
    if(p[cb(a,b,0)]>=max(p[cb(a,0,c)],p[cb(0,b,c)])){
      zmin[cb1(a,b,c)]<-ztest(C@n[cb(a,b,c)],C@n[cb(a,b,0)],p[cb(a,b,c)],p[cb(a,b,0)])
    }
    if(p[cb(a,0,c)]>=max(p[cb(a,b,0)],p[cb(0,b,c)])){
      zmin[cb1(a,b,c)]<-ztest(C@n[cb(a,b,c)],C@n[cb(a,0,c)],p[cb(a,b,c)],p[cb(a,0,c)])
    }
    if(p[cb(0,b,c)]>=max(p[cb(a,b,0)],p[cb(a,0,c)])){
      zmin[cb1(a,b,c)]<-ztest(C@n[cb(a,b,c)],C@n[cb(0,b,c)],p[cb(a,b,c)],p[cb(0,b,c)])
    }
  }}}
  zave<-mean(zmin)/sqrt(var(zmin)/(C@D[1]*C@D[2]*C@D[3]))
  if(is.null(simerror)) simerror<-9 #This is necessary for interfacing to C++
  if(is.null(nboot)) nboot<-900
  results<-.Call("avebinomial3",C@n,p,c(nboot,C@D),c(zave,simerror),PACKAGE="bifactorial")
  simerror<-sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3))
  dauer<-proc.time()[3]-dauer
  new("avetest",stat=zave,p=round(results[[1]]/results[[2]],4),
      nboot=results[[2]],simerror=simerror,method="Bootstrap",duration=dauer)
}
#MAX-test based on a Z-statistic in trifactorial designs
maxbinomial3Boot<-function(C,nboot=NULL,simerror=NULL,...){
  pv<-mintest(C,test="ztest",nboot=nboot,simerror=simerror,...)
  pmax<-round(min(as.numeric(pv@p)),4)
  new("maxtest",p=pmax,stat=max(pv@stat),name=pv@gnames[pv@p==pmax],nboot=pv@nboot,simerror=pv@simerror,method="bootstrap",duration=pv@duration)
}
