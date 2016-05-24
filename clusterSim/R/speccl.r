.GDMKernel <-function (gdm,sigma = 1) 
{
  resul<-exp(as.matrix(-sigma*gdm))  
  for(i in 1:ncol(gdm)){
    resul[i,i]<-0
  }
  resul
}



.GausKernel<-function(d,sigma=1){
  resul<-exp(as.matrix(-sigma*d))  
  for(i in 1:ncol(resul)){
    resul[i,i]<-0
  }
  resul
}

.ddist<-function(dane,distType){
  res<-NULL
  if(distType=="GDM1"){
    res<-GDM1(dane)
  }
  else if(distType=="GDM2"){
    res<-GDM2(dane)
  }
  else if(distType=="sEuclidean"){
    res<-dist(dane)^2
  }
  else if(distType=="BC"){
    res<-dist.BC(dane)
  }
  else if(distType=="SM"){
    res<-dist.SM(dane)
  }
  else{
    res<-dist(dane,method=distType)
  }
  res
}


speccl<-function(data,nc,distance="GDM1",sigma="automatic",sigma.interval="default",mod.sample=0.75,R=10,iterations=3,na.action=na.omit,...)
{
  if(sigma=='automatic'){
    sigmaSimulation<-TRUE
  }
  else{
    sigmaSimulation<-FALSE
    sigma<-as.numeric(sigma)
  }
  DEBUG<-FALSE
  globalOk<-FALSE
  silDebug=TRUE
  badSigma<-NULL
  tries<-0
  while(!globalOk){
  step<-0
  sigWithinss<--1
  ok<-FALSE
  while(!ok && step<6){
  step<-step+1
  #print(paste("step",step))
  x<-data
  bootstrap<-x[sample(1:nrow(x),nrow(x)*mod.sample),]
  levelsPower<-R;
    levels<-iterations
    lstart<-0
    lend<-sum(.ddist(x,distance)) # tutaj suma odleglosci (euklidesowa, miejska, kwadrat euklidesowej)
    if(distance=="sEuclidean"){ lend<-sqrt(lend) }# inna gorna granica sigmy , w tym przypadku jest pierwiastek
    if(sigma.interval!="default") {
      lend<-sigma.interval
     }
    # zeby to byla suma odl. euklidesowej podstawienie :   lend<-.ddist(x,"euclidean")
  lby<-lend
  lstartend<-lend
  sig<-sample(1:lstartend,1)
  if(sigmaSimulation){
    for(ll in levels:1){
      lby<-lby/levelsPower
      sigmas<-(seq(lstart,lend-lby,by=lby)+seq(lstart+lby,lend,by=lby))/2
      oldsigma<-sig
      i<-0
      for (sigma in sigmas){
         #print(paste("iteration level=",ll,"sigma=",sigma))
        if(distance=="GDM1" || distance=="GDM2"){
          ka<-.GDMKernel(as.matrix(dist.GDM(bootstrap,method=distance)),sigma)
        }
        else if(distance=="BC"){
          ka<-.GausKernel(as.matrix(dist.BC(bootstrap)),sigma)
        }
        else if(distance=="SM"){
          ka<-.GausKernel(as.matrix(dist.SM(bootstrap)),sigma)
        }
        else if(distance=="sEuclidean"){
          ka<-.GausKernel(as.matrix(dist(bootstrap))^2,sigma)
        }
        else{
          dd<-try(dist(bootstrap,method=distance),silent=silDebug)
          if(class(dd)=="try-error"){
            dd<-try(dist.binary(bootstrap,method=distance),silent=silDebug)
          }
          if(class(dd)=="try-error"){
            stop(paste("unknown distance method ",distance))
          }
          ka<-.GausKernel(as.matrix(dd),sigma)
        }
        d<-1/sqrt(rowSums(ka))
        l<-d * ka %*% diag(d)
        ei<-NULL
        tf<-function(l,nc){eigen(l,symmetric=TRUE)$vectors[,1:nc]}
        ei<-try(tf(l,nc),silent=silDebug)
        #bbootstrap<<-bootstrap
        #dd<<-d
        #kka<<-ka
        #ll<<-l
        #ssigma<<-sigma
        ##print(class(ei))
        if(class(ei)!="try-error"){
          if(!is.null(ei)  && is.numeric(ei)){
            yi<-try(ei/sqrt(rowSums(ei^2)),silent=silDebug)
            if(sum(is.na(yi))==0){
            res<-try(kmeans(yi, yi[initial.Centers(yi,nc),],...),silent=silDebug)
            if(class(res)=="try-error"){
              res<-list(withinss=1e10)
              next
            }
            if(sum(res$withinss)<sigWithinss || sigWithinss==-1){
              ok<-TRUE
              sig<-sigma
              sigWithinss<-sum(res$withinss)
            }
          }
          i<-i+1
          }
          else{
            na.action(ei)
          }
        }
        else{
            #stop("BAD EIGEN")
        }
      }
      #ssig<<-sig
      #ooldsigma<<-oldsigma
      if(is.null(sig) || (!is.null(oldsigma) && oldsigma==sig)){
        lstart<-lstart/R
        lend<-lend/R
      }
      else{
        lstart<-sig-0.5*lby
        lend<-sig+0.5*lby
      }
    }
  }
  else{
    ok<-TRUE
    sig<-as.numeric(sigma)
  }
  }
  #print(paste("sigma found",sig))
  if(step>=6){
    sig<-sample(1:lstartend,1)
    if(distance=="manhattan")    sig<-sample(1:10,1)
    #print("step six or higer")
  }
  if(!is.null(badSigma)){
    for(ss in badSigma){
      if(abs(sig-ss)<0.5){
        sig<-sample(1:lstartend,1)
        if(distance=="manhattan")    sig<-sample(1:10,1)
        #print(paste("new random sigma",sig))
      }
    }
  }
  
  globalOk<-TRUE
  if(distance=="GDM1" || distance=="GDM2"){
    scdist<-dist.GDM(x,method=distance)
    km<-.GDMKernel(as.matrix(scdist),sig)
  }
  else if(distance=="BC"){
    scdist<-dist.BC(x)
    km<-.GausKernel(as.matrix(scdist),sig)
  }
  else if(distance=="SM"){
    scdist<-dist.SM(x)
    km<-.GausKernel(as.matrix(scdist),sig)
  }
  else if(distance=="sEuclidean"){
    scdist<-dist(x)^2
    km<-.GausKernel(as.matrix(scdist),sig)
  }
  else{
        scdist<-try(dist(x,method=distance),silent=silDebug)
        if(class(dd)=="try-error"){
          scdist<-try(dist.binary(x,method=distance),silent=silDebug)
        }
        if(class(scdist)=="try-error"){
          stop(paste("unknown distance method ",distance))
          globalOk<-FALSE
        }
        km<-.GausKernel(as.matrix(scdist),sig)
  }
  
  diag(km)<-0
  d<-1/sqrt(rowSums(km))
  l<-d * km %*% diag(d)
  if(getRversion() >= '3.0'){
  ei<-try(eigen(l,symmetric=T)$vectors[, 1:nc],silent=silDebug)
  }
  else{
  ei<-try(eigen(l)$vectors[, 1:nc],silent=silDebug)
  }
  if(class(ei)=="try-error"){
    #stop(paste("Not possible to calculate eigenvalues, try with other distance type - ",distance))
    globalOk<-FALSE
    #print("bad eigen")
  }
  if(globalOk){yi<-ei/sqrt(rowSums(ei^2))}
  if(globalOk){
    res<-try(kmeans(yi, yi[initial.Centers(yi, nc),], ...),silent=silDebug)
    if(class(res)=="try-error"){
      res<-try(kmeans(yi,nc,...),silent=silDebug)
    }
  }
  if(globalOk && class(res)=="try-error"){
    #yyi<<-yi
    #print("bad clustering")
    if(is.character(all.equal(na.action,na.omit))){     # if all.equals not returns true it returns string"
      tries<-tries+1
      if(tries<5){
        stop(paste("Not possible to do clustering, try with other distance type - ",distance))
      }
    }
    globalOk<-FALSE
  }
  if(!globalOk){
    badSigma<-c(badSigma,sig)}
    #print(paste("Bad sigma",badSigma))
  }
  return(list(clusters = res$cluster, size = res$size,withinss=res$withins,sigma=sig,Ematrix=ei,Ymatrix=yi,scdist=scdist))
}

#data(data_binary)
##print(speccl(data_binary,nc=3,distance=1,sigma="automatic",mod.sample=0.75,R=10,iterations=3))
#grnd2<-cluster.Gen(50,model=4,dataType="m",numNoisyVar=1)
#data<-as.matrix(grnd2$data)
#colornames<-c("red","blue","green")
#grnd2$clusters[grnd2$clusters==0]<-length(colornames)
#plot(grnd2$data,col=colornames[grnd2$clusters])
#res2<-speccl(data,nc=3,distance="manhattan",sigma="automatic",sigma.interval="default",mod.sample=0.75,R=10,iterations=3)
##print(res2$sigma)
#cRand<-comparing.Partitions(grnd2$clusters,res2$clusters,type="crand")
##print(cRand)

