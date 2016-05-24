# #######################################################################################
# ### sam stuff for Kobe ################################################################
# #######################################################################################

utils::globalVariables(c("pctl","FLQuant","stock.n<-","propagate","stock.n","harvest<-","mvrnorm"))

readSam<-function(file, reduced=FALSE){
  # Function to read a basic fit
  ret<-list()
  parfile<-as.numeric(scan(paste(file,'.par', sep=''), 
                           what='', n=16, quiet=TRUE)[c(6,11,16)])
  ret$nopar<-as.integer(parfile[1])
  ret$nlogl<-parfile[2]
  ret$maxgrad<-parfile[3]
  rep<-scan(paste(file,'.rep', sep=''), quiet=TRUE)
  ret$res<-read.table(paste(file,'.res', sep=''),header=FALSE)
  ret$stateDim<-rep[1]
  ret$years<-rep[-1]
  file<-paste(file,'.cor', sep='')
  lin<-readLines(file)
  ret$npar<-length(lin)-2
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2])
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!=''])
  ret$names<-unlist(lapply(sublin,function(x)x[2]))
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3])))
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4])))
  
  ret$cor<-matrix(NA, ret$npar, ret$npar)
  for(i in 1:ret$npar){
    ret$cor[1:i,i]<-as.numeric(unlist(lapply(sublin[i],
                                             function(x)x[5:(4+i)])))
    ret$cor[i,1:i]<-ret$cor[1:i,i]
  }
  ret$cov<-ret$cor*(ret$std%o%ret$std)
  mslh<-function(name){
    idx<-which(ret$names==name)
    x<-cbind(ret$est[idx], ret$std[idx], ret$est[idx]-2*ret$std[idx], 
             ret$est[idx]+2*ret$std[idx])
    colnames(x)<-c('est', 'std', 'low', 'hig')
    return(x)
  }
  ret$ssb<-mslh('ssb')
  ret$fbar<-mslh('fbar')
  ret$tsb<-mslh('tsb')
  ret$logssb<-mslh('logssb')
  ret$logfbar<-mslh('logfbar')
  ret$logtsb<-mslh('logtsb')
  ret$logscale<-mslh('logScale')
  ret$logFpar<-mslh('logFpar')
  ret$logCatch<-mslh('logCatch')
  x<-mslh('U')
  ret$stateEst<-matrix(x[,1],ncol=ret$stateDim, byrow=TRUE)
  ret$stateStd<-matrix(x[,2],ncol=ret$stateDim, byrow=TRUE)
  ret$stateLow<-matrix(x[,3],ncol=ret$stateDim, byrow=TRUE)
  ret$stateHig<-matrix(x[,4],ncol=ret$stateDim, byrow=TRUE)
  ret$R<-cbind(exp(ret$stateEst[,1]), NA, exp(ret$stateLow[,1]), 
               exp(ret$stateHig[,1]))
  
  if(reduced)
    ret<-ret[which(!names(ret)%in%c('cov','cor'))]
  
  return(ret)}

kobeSamFn=function(object,what=c("trks","smry","pts","ellipse")[1],prob=c(0.75,0.5,.25),nits=1000,bmsy=1,fmsy=1){

    #require(ellipse)
  
    res   =readSam(object)
    
    trks=cbind(melt(data.frame(year=res$years,res$ssb),id="year"),
               melt(res$fbar))[,-c(4,5)]
    names(trks)[2:4]=c("quantity","stock","harvest")
    
    posSsb =seq(length(res$names))[res$names=="ssb"]
    posFbar=seq(length(res$names))[res$names=="fbar"]
    
    t.=trks[trks$quantity=="est",]
    t.=t.[dim(t.)[1],]
    
    currentCov=res$cov[c(rev(posFbar)[1],rev(posSsb)[1]),
                       c(rev(posFbar)[1],rev(posSsb)[1])]
    ellipse=as.data.frame(ellipse(currentCov,level=prob[1]))
    
    names(ellipse)=c("harvest","stock")
    ellipse=transform(ellipse, stock=stock+t.$stock,harvest=harvest+t.$harvest)
    
    rnd=c(t(mvrnorm(nits,res$est[c(rev(posFbar)[1],rev(posSsb)[1])],currentCov)))
    rnd=as.data.frame(t(array(c(rnd),c(2,nits))))
    names(rnd)=c("harvest","stock")
 
    res=llply(list(trks=trks,ellipse=ellipse,pts=rnd), function(x,bmsy,fmsy) transform(x, stock=stock/bmsy, harvest=harvest/fmsy), bmsy=bmsy, fmsy=fmsy)
    
    return(res)
    
}

#object="/home/laurie/Desktop/ICCAT/SCRS/2013/SWO/SAM/swon/run/sam"
#t.=kobeSam(object)
#kobePhase(t.$ellipse)+geom_path(aes(stock,harvest))+geom_path(aes(stock,harvest),data=subset(t.$trks,quantity=="est"))+geom_point(aes(stock,harvest),data=t.$pts)
 
#### exported function
kobeSam=function(object,what=c("trks","smry","pts","ellipse")[1],prob=c(0.75,0.5,.25),nits=1000,bmsy=1,fmsy=1){
    
  if (length(object)>1){
     res=kobeSamFn(object,what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy)
  } else {  
     res=mlply(object,   function(x,what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy)
                        kobeSamFn(x,what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy),
                                    what=what,prob=prob,nits=nits,bmsy=bmsy,fmsy=fmsy)
                    
      res=list(trks   =ldply(res, function(x) x$trks),
               pts    =ldply(res, function(x) x$pts),
               smry   =ldply(res, function(x) x$smry),
               ellipse=ldply(res, function(x) x$ellipse))
      }
    
    if (length(what)==1) return(res[[what]]) else return(res[what]) 
    }

samNF=function(sam,dir,nits=1000){
  res=readSam(dir)
  
  yrs=as.numeric(res$years)
  ags=min(res$res[,3]):max(res$res[,3])
  mm= res$names=="U"
  
  rnd=c(t(mvrnorm(nits,res$est[mm],res$cov[mm,mm])))
  rnd=array(c(rnd),c(length(ags)*2-1,length(yrs),nits))
  rnd=FLQuant(exp(c(rnd)),dimnames=list(age=seq(length(ags)*2-1),year=yrs,iter=seq(nits)))
   
  stock.n(sam)=propagate(stock.n(sam),nits)
  harvest(sam)=propagate(harvest(sam),nits)
   
   stock.n(sam)[                  ]=rnd[seq(ags)]
   harvest(sam)[seq(length(ags)-1)]=rnd[seq(length(ags)-1)+length(ags)]
   harvest(sam)[max(ags)          ]=rnd[length(ags)*2-1]
   
   return(sam)}
