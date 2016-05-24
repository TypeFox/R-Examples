maximinSLHD <-
function(t,m,k,power=15,nstarts=1,itermax=100,total_iter=1000000){
  
  n=m*t
  
  t00<-Sys.time()
  aaa<-.C("maximinSLHD",as.integer(m),as.integer(k),as.integer(t),as.integer(power),as.integer(nstarts),
          as.integer(itermax),as.integer(total_iter),design=integer(m*t*k),measure=double(1), 
          temp0=double(1),PACKAGE="SLHD")
  t01<-Sys.time()
  
  time_rec=t01-t00
  
  dd=matrix(aaa$design,ncol=k,nrow=n,byrow=TRUE)
  dds<-(apply(dd,2,rank)-0.5)/n
    
  if(t>1){
    Slice<-rep(1:t,each=m)
    dd<-cbind(Slice,dd)
    dds<-cbind(Slice,dds)
  }
  
  val<-list(Design=dd,measure=aaa$measure,StandDesign=dds,temp0=aaa$temp0,time_rec=time_rec)
  
  return(val)
}
