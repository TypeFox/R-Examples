#Implementations of bootstrap algorithms for the min-test (Laska and Meisner, 1989)
#min-test based on Student's t-test in bifactorial designs
pvstudent2Boot<-function(C,nboot,simerror,...){
  if(is.binary(C@data)) warning("For binary data, the Z-statistic should be used.")
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  dauer<-proc.time()[3] #Measure duration of simulation
  Y<-numeric(0)   
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){Y<-c(Y,C@data[[cb(a,b)]])}}
  tmin<-matrix(nrow=C@D[1],ncol=C@D[2])
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){
    if(mean(C@data[[cb(a,0)]])>=mean(C@data[[cb(0,b)]])){
      tmin[a,b]<-ttest(C@data[[cb(a,b)]],C@data[[cb(a,0)]])
    }
    else{tmin[a,b]<-ttest(C@data[[cb(a,b)]],C@data[[cb(0,b)]])}
  }}
  results<-.Call("student2",Y,matrix(C@n,nrow=C@D[1]+1,byrow=TRUE),tmin,
                 c(nboot,C@D[1],C@D[2]),simerror,PACKAGE="bifactorial")
  simerror<-sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3))
  dauer<-proc.time()[3]-dauer
  return(new("mintest",gnames=kombi(C@D),method="Bootstrap",
             p=round(results[[1]]/results[[2]],4),stat=as.numeric(t(tmin)),
             nboot=results[[2]],simerror=simerror,duration=dauer))
}
#min-test based on Student's t-test in trifactorial designs
pvstudent3Boot<-function(C,nboot,simerror,...){
  if(is.binary(C@data)) warning("For binary data, the Z-statistic should be used.")
  cb<-function(a,b,c){(a*(C@D[2]+1)*(C@D[3]+1))+(b*(C@D[3]+1))+c+1}
  cb1<-function(a,b,c){((a-1)*(C@D[2])*(C@D[3]))+((b-1)*(C@D[3]))+c}
  dauer<-proc.time()[3] #Measure duration of simulation
  tmin<-Y<-numeric(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){for(c in 0:C@D[3]){Y<-c(Y,C@data[[cb(a,b,c)]])}}}
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){for(c in 1:C@D[3]){
    if(mean(C@data[[cb(0,b,c)]])>=max(mean(C@data[[cb(a,0,c)]]),mean(C@data[[cb(a,b,0)]]))){
      tmin[cb1(a,b,c)]<-ttest(C@data[[cb(a,b,c)]],C@data[[cb(0,b,c)]])
    }
    if(mean(C@data[[cb(a,0,c)]])>=max(mean(C@data[[cb(0,b,c)]]),mean(C@data[[cb(a,b,0)]]))){
      tmin[cb1(a,b,c)]<-ttest(C@data[[cb(a,b,c)]],C@data[[cb(a,0,c)]])
    }
    if(mean(C@data[[cb(a,b,0)]])>=max(mean(C@data[[cb(a,0,c)]]),mean(C@data[[cb(0,b,c)]]))){
      tmin[cb1(a,b,c)]<-ttest(C@data[[cb(a,b,c)]],C@data[[cb(a,b,0)]])
    }
  }}}
  results<-.Call("student3",Y,C@n,tmin,c(nboot,C@D),simerror,PACKAGE="bifactorial")
  simerror<-sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3))
  dauer<-proc.time()[3]-dauer
  return(new("mintest",gnames=kombi(C@D),method="Bootstrap",
             p=round(results[[1]]/results[[2]],4),stat=as.numeric(t(tmin)),
             nboot=results[[2]],simerror=simerror,duration=dauer))
}
#min-test based on a Z-statistic in bifactorial designs
pvbinomial2Boot<-function(C,nboot,simerror,...){
  if(!is.binary(C@data)) stop("This test needs binary data.")
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  dauer<-proc.time()[3] #Measure duration of simulation
  n<-matrix(C@n,nrow=C@D[1]+1,byrow=TRUE)
  p<-matrix(nrow=C@D[1]+1,ncol=C@D[2]+1)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){p[a+1,b+1]=sum(C@data[[cb(a,b)]])/n[a+1,b+1]}}
  zmin<-matrix(nrow=C@D[1],ncol=C@D[2])
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){
    zmin[a,b]<-if(p[a+1,1]>=p[1,b+1]){ztest(n[a+1,b+1],n[a+1,1],p[a+1,b+1],p[a+1,1])}
    else{ztest(n[a+1,b+1],n[1,b+1],p[a+1,b+1],p[1,b+1])}
  }}
  results<-.Call("binomial2",n,p,zmin,c(nboot,C@D),simerror,PACKAGE="bifactorial")
  simerror<-sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3))
  dauer<-proc.time()[3]-dauer
  return(new("mintest",gnames=kombi(C@D),method="Bootstrap",
             p=round(results[[1]]/results[[2]],4),stat=as.numeric(t(zmin)),
             nboot=results[[2]],simerror=simerror,duration=dauer))
}
#min-test based on a Z-statistic in trifactorial designs
pvbinomial3Boot<-function(C,nboot,simerror,...){
  if(!is.binary(C@data)) stop("This test needs binary data.")
  cb<-function(a,b,c){(a*(C@D[2]+1)*(C@D[3]+1))+(b*(C@D[3]+1))+c+1}
  cb1<-function(a,b,c){((a-1)*(C@D[2])*(C@D[3]))+((b-1)*(C@D[3]))+c}
  dauer<-proc.time()[3] #Measure duration of simulation
  p<-zmin<-numeric(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){for(c in 0:C@D[3]){
    p[cb(a,b,c)]=mean(C@data[[cb(a,b,c)]])
  }}}
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){for(c in 1:C@D[3]){
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
  results<-.Call("binomial3",C@n,p,zmin,c(nboot,C@D),simerror,PACKAGE="bifactorial")
  simerror<-sqrt((results[[1]]/results[[2]]^2)-(results[[1]]^2/results[[2]]^3))
  dauer<-proc.time()[3]-dauer
  return(new("mintest",gnames=kombi(C@D),method="Bootstrap",
             p=round(results[[1]]/results[[2]],4),stat=as.numeric(t(zmin)),
             nboot=results[[2]],simerror=simerror,duration=dauer))
}
