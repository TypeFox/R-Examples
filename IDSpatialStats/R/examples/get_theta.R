data(DengueSimulationR02)

r.max<-seq(20,1000,20)
r.min<-seq(0,980,20)

sero.type.func<-function(a,b,tlimit=20){
  if(a[5]==b[5]&(abs(a[3]-b[3])<=tlimit)){rc=1}
  else{rc=2}
  return(rc)
}

sero.theta<-get.theta(DengueSimR02,sero.type.func,r=r.max,r.low=r.min)
