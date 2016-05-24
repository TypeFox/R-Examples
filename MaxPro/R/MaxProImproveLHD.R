 

MaxProImproveLHD<-
function(InitialDesign,s=2,localopm=FALSE,temp0=0,nstarts=1,itermax=400,total_iter=1000000){
    
  DD<-as.matrix(InitialDesign)
  #DD<-apply(DD,2,rank)
  m<-nrow(DD)
  k<-ncol(DD)
  n=m
  
  if(temp0==0){
    width<-mean(apply(DD,2,range))
    avgdist1=1/(width)
    avgdist2=(1/(width^(k-1)*(width-width/(n-1))))^(1/k)
    delta=avgdist2-avgdist1
    temp0=-delta/log(0.99)
  }
  
  localopm<-as.numeric(localopm)
 
  predesign=c(t(DD))
    
  t00<-Sys.time()
  aaa<-.C("MaxProImprove",as.integer(m),as.integer(k),as.integer(localopm),as.double(predesign),as.integer(nstarts),
          as.integer(itermax),as.integer(total_iter),design=double(m*k),measure=double(1), 
          as.double(temp0),ntotalI=integer(1),as.integer(s),PACKAGE="MaxPro")
  t01<-Sys.time()
  
  time_rec=t01-t00  
  Design=matrix(aaa$design,ncol=k,nrow=m,byrow=TRUE)
      
  val<-list(Design=Design,temp0=temp0,measure=aaa$measure,time_rec=time_rec,ntotal=aaa$ntotalI)
  
  return(val)
}
