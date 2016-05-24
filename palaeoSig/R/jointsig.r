jointsig<-function(spp,fos,var1, var2,method="randomTF",n=99,r=32, ...){
  if(r%%4!=0) warning("r not divisible by 4")
  v1<-sin(seq(0,2*pi,length=r+1))[-(r+1)]
  v2<-cos(seq(0,2*pi,length=r+1))[-(r+1)]
  syn.env<-as.data.frame(mapply(function(v1,v2){scale(var1)*v1+scale(var2)*v2},v1=v1,v2=v2))
  if(method=="randomTF"){
    res<-randomTF(spp,fos,..., env=syn.env, n=n)
  }else if(method=="obs.cor"){
    warning("obscor can pathological results in jointsig")
    res<-list()
    
    res$EX<-sapply(syn.env,function(e)obs.cor(spp,fos,..., env=e, n=0)$ob$res$wc) ######edit to match new obscor
    res$sim.ex<-obs.cor(spp,fos,...,env=syn.env[,1], n=n)$sim[,2]
 
  }
  res$v1<-v1
  res$v2<-v2
  class(res)<-"js"
  return(res)
}



