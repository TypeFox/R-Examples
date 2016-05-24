#Implementations of bootstrap algorithms for simulataneous confidence intervals
#Metric data: Intervals based on Student's t-test in bifactorial designs
intstudent2Boot<-function(C,nboot,simerror,alpha,...){
  dauer<-proc.time()[3]
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  cb1<-function(a,b){((a-1)*C@D[2])+b}
  kiu<-kio<-Y<-mx<-vx<-numeric(0)
  cnames<-character(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){
    mx[cb(a,b)]=mean(C@data[[cb(a,b)]])
    vx[cb(a,b)]=var(C@data[[cb(a,b)]])
  }}
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){
   Y<-c(Y,C@data[[cb(a,b)]]-mx[cb(a,b)])
  }}
  if(is.null(nboot)) nboot<-1
  if(!is.null(simerror)) nboot=max(nboot,alpha*(1-alpha)/(simerror^2))
  results<-.Call("kritstudent2",Y,C@n,c(nboot,C@D),PACKAGE="bifactorial")
  qu<--quantile(results[[2]],1-alpha/2)
  qo<-quantile(results[[2]],1-alpha/2)
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){
    kiu<-c(kiu,mx[cb(a,b)]-mx[cb(a,0)]+qu*sqrt((vx[cb(a,b)]/C@n[cb(a,b)])+(vx[cb(a,0)]/C@n[cb(a,0)])))
    kiu<-c(kiu,mx[cb(a,b)]-mx[cb(0,b)]+qu*sqrt((vx[cb(a,b)]/C@n[cb(a,b)])+(vx[cb(0,b)]/C@n[cb(0,b)])))
    kio<-c(kio,mx[cb(a,b)]-mx[cb(a,0)]+qo*sqrt((vx[cb(a,b)]/C@n[cb(a,b)])+(vx[cb(a,0)]/C@n[cb(a,0)])))
    kio<-c(kio,mx[cb(a,b)]-mx[cb(0,b)]+qo*sqrt((vx[cb(a,b)]/C@n[cb(a,b)])+(vx[cb(0,b)]/C@n[cb(0,b)])))
    cnames<-c(cnames,paste("(",a,",",b,")-(",a,",0)",sep=""),paste("(",a,",",b,")-(0,",b,")",sep=""))
  }}
  dauer=proc.time()[3]-dauer
  new("margint",cnames=cnames,kiu=kiu,kio=kio,method="Bootstrap",nboot=nboot,
      simerror=sqrt(alpha*(1-alpha)/nboot),alpha=alpha,duration=dauer)
}
#Intervals based on Student's t-test in trifactorial designs
intstudent3Boot<-function(C,nboot,simerror,alpha,...){
  dauer<-proc.time()[3]
  cb<-function(a,b,c){(a*(C@D[2]+1)*(C@D[3]+1))+(b*(C@D[3]+1))+c+1}
  cb1<-function(a,b,c){((a-1)*(C@D[2])*(C@D[3]))+((b-1)*(C@D[3]))+c}
  kiu<-kio<-Y<-mx<-vx<-numeric(0)
  cnames<-character(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){for(c in 0:C@D[3]){
    mx[cb(a,b,c)]=mean(C@data[[cb(a,b,c)]])
    vx[cb(a,b,c)]=var(C@data[[cb(a,b,c)]])
  }}}
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){for(c in 0:C@D[3]){
    Y<-c(Y,C@data[[cb(a,b,c)]]-mean(C@data[[cb(a,b,c)]]))
  }}}
  if(is.null(nboot)) nboot<-1
  if(!is.null(simerror)) nboot=max(nboot,alpha*(1-alpha)/(simerror^2))
  results<-.Call("kritstudent2",Y,C@n,c(nboot,C@D),PACKAGE="bifactorial")
  qu<--quantile(results[[2]],1-alpha/2)
  qo<-quantile(results[[2]],1-alpha/2)
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){for(c in 1:C@D[3]){
    kiu<-c(kiu,mx[cb(a,b,c)]-mx[cb(a,b,0)]+qu*sqrt((vx[cb(a,b,c)]/C@n[cb(a,b,c)])+(vx[cb(a,b,0)]/C@n[cb(a,b,0)])))
    kiu<-c(kiu,mx[cb(a,b,c)]-mx[cb(a,0,c)]+qu*sqrt((vx[cb(a,b,c)]/C@n[cb(a,b,c)])+(vx[cb(a,0,c)]/C@n[cb(a,0,c)])))
    kiu<-c(kiu,mx[cb(a,b,c)]-mx[cb(0,b,c)]+qu*sqrt((vx[cb(a,b,c)]/C@n[cb(a,b,c)])+(vx[cb(0,b,c)]/C@n[cb(0,b,c)])))
    kio<-c(kio,mx[cb(a,b,c)]-mx[cb(a,b,0)]+qo*sqrt((vx[cb(a,b,c)]/C@n[cb(a,b,c)])+(vx[cb(a,b,0)]/C@n[cb(a,b,0)])))
    kio<-c(kio,mx[cb(a,b,c)]-mx[cb(a,0,c)]+qo*sqrt((vx[cb(a,b,c)]/C@n[cb(a,b,c)])+(vx[cb(a,0,c)]/C@n[cb(a,0,c)])))
    kio<-c(kio,mx[cb(a,b,c)]-mx[cb(0,b,c)]+qo*sqrt((vx[cb(a,b,c)]/C@n[cb(a,b,c)])+(vx[cb(0,b,c)]/C@n[cb(0,b,c)])))
    cnames<-c(cnames,paste("(",a,",",b,",",c,")-(",a,",",b,",0)",sep=""),
              paste("(",a,",",b,",",c,")-(",a,",0,",c,")",sep=""),paste("(",a,",",b,",",c,")-(0,",b,",",c,")",sep=""))
  }}}
  dauer=proc.time()[3]-dauer
  new("margint",cnames=cnames,kiu=kiu,kio=kio,method="Bootstrap",nboot=nboot,
      simerror=sqrt(alpha*(1-alpha)/nboot),alpha=alpha,duration=dauer)
}
#Binary data: Intervals based on a Z-statistic in bifactorial designs
intbinomial2Boot<-function(C,nboot,simerror,alpha,...){
  cb<-function(a,b){(a*(C@D[2]+1))+b+1}
  cb1<-function(a,b){((a-1)*C@D[2])+b}
  v<-function(x) x*(1-x)
  dauer<-proc.time()[3]
  p<-kiu<-kio<-mx<-Y<-vx<-numeric(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){
    p[cb(a,b)]=mean(C@data[[cb(a,b)]])
  }}
  if(is.null(nboot)) nboot<-1
  if(!is.null(simerror)) nboot=max(nboot,round(alpha*(1-alpha)/(simerror^2)))
  results<-.Call("kritbinomial2",p,C@n,c(nboot,C@D),PACKAGE="bifactorial")
  qu<--quantile(results[[2]],1-alpha/2)
  qo<-quantile(results[[2]],1-alpha/2)
  cnames<-character(0)
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){
    kiu<-c(kiu,p[cb(a,b)]-p[cb(a,0)]+qu*sqrt((v(p[cb(a,b)])/C@n[cb(a,b)])+(v(p[cb(a,0)])/C@n[cb(a,0)])))
    kiu<-c(kiu,p[cb(a,b)]-p[cb(0,b)]+qu*sqrt((v(p[cb(a,b)])/C@n[cb(a,b)])+(v(p[cb(0,b)])/C@n[cb(0,b)])))
    kio<-c(kio,p[cb(a,b)]-p[cb(a,0)]+qo*sqrt((v(p[cb(a,b)])/C@n[cb(a,b)])+(v(p[cb(a,0)])/C@n[cb(a,0)])))
    kio<-c(kio,p[cb(a,b)]-p[cb(0,b)]+qo*sqrt((v(p[cb(a,b)])/C@n[cb(a,b)])+(v(p[cb(0,b)])/C@n[cb(0,b)])))
    cnames<-c(cnames,paste("(",a,",",b,")-(",a,",0)",sep=""),paste("(",a,",",b,")-(0,",b,")",sep=""))
  }}
  dauer=proc.time()[3]-dauer
  new("margint",cnames=cnames,kiu=kiu,kio=kio,method="Bootstrap",nboot=nboot,
      simerror=sqrt(alpha*(1-alpha)/nboot),alpha=alpha,duration=dauer)
}
#Intervals based on a Z-statistic in trifactorial designs
intbinomial3Boot<-function(C,nboot,simerror,alpha,...){
  cb<-function(a,b,c){(a*(C@D[2]+1)*(C@D[3]+1))+(b*(C@D[3]+1))+c+1}
  cb1<-function(a,b,c){((a-1)*(C@D[2])*(C@D[3]))+((b-1)*(C@D[3]))+c}
  v<-function(x) x*(1-x)
  dauer<-proc.time()[3] 
  p<-kiu<-kio<-mx<-Y<-vx<-numeric(0)
  for(a in 0:C@D[1]){for(b in 0:C@D[2]){for(c in 0:C@D[3]){
    p[cb(a,b,c)]=mean(C@data[[cb(a,b,c)]])
  }}}
  if(is.null(nboot)) nboot<-1
  if(!is.null(simerror)) nboot=max(nboot,round(alpha*(1-alpha)/(simerror^2)))
  results<-.Call("kritbinomial3",p,C@n,c(nboot,C@D),PACKAGE="bifactorial")
  qu<-1-quantile(results[[2]],1-alpha/2)
  qo<-quantile(results[[2]],1-alpha/2)
  cnames<-character(0)
  for(a in 1:C@D[1]){for(b in 1:C@D[2]){for(c in 1:C@D[3]){
    kiu<-c(kiu,p[cb(a,b,c)]-p[cb(a,b,0)]+qu*sqrt((v(p[cb(a,b,c)])/C@n[cb(a,b,c)])+(v(p[cb(a,b,0)])/C@n[cb(a,b,0)])))
    kiu<-c(kiu,p[cb(a,b,c)]-p[cb(a,0,c)]+qu*sqrt((v(p[cb(a,b,c)])/C@n[cb(a,b,c)])+(v(p[cb(a,0,c)])/C@n[cb(a,0,c)])))
    kiu<-c(kiu,p[cb(a,b,c)]-p[cb(0,b,c)]+qu*sqrt((v(p[cb(a,b,c)])/C@n[cb(a,b,c)])+(v(p[cb(0,b,c)])/C@n[cb(0,b,c)])))
    kio<-c(kio,p[cb(a,b,c)]-p[cb(a,b,0)]+qo*sqrt((v(p[cb(a,b,c)])/C@n[cb(a,b,c)])+(v(p[cb(a,b,0)])/C@n[cb(a,b,0)])))
    kio<-c(kio,p[cb(a,b,c)]-p[cb(a,0,c)]+qo*sqrt((v(p[cb(a,b,c)])/C@n[cb(a,b,c)])+(v(p[cb(a,0,c)])/C@n[cb(a,0,c)])))
    kio<-c(kio,p[cb(a,b,c)]-p[cb(0,b,c)]+qo*sqrt((v(p[cb(a,b,c)])/C@n[cb(a,b,c)])+(v(p[cb(0,b,c)])/C@n[cb(0,b,c)])))
    cnames<-c(cnames,paste("(",a,",",b,",",c,")-(",a,",",b,",0)",sep=""),
              paste("(",a,",",b,",",c,")-(",a,",0,",c,")",sep=""),paste("(",a,",",b,",",c,")-(0,",b,",",c,")",sep=""))
  }}}
  dauer=proc.time()[3]-dauer
  new("margint",cnames=cnames,kiu=kiu,kio=kio,method="Bootstrap",nboot=nboot,
      simerror=sqrt(alpha*(1-alpha)/nboot),alpha=alpha,duration=dauer)
}
