classifyEEG <-
function(y,data){
  
  LL<-y$L
  nel=y$nch
  NN<-nrow(data)
  if(NN<LL) data<-rbind(matrix(rep(0,(LL-NN)*nel),ncol=nel),data) else {
    data<-data[(NN-LL+1):NN,]
  }
  featype<-y$featype


  
  #####SPECTRO
  if (sum(c("f3","f4")%in%featype)>0) {
    specdata<-.spec.pgram(data)
  }

  ####PCA
  if (sum(c("f2","f5")%in%featype)>0) {
    ncomps<-y$ncomps
    rot<-eigen(var(data))$vectors[,1:ncomps]
    loaddata<-rot
    pcadata<-as.matrix(data)%*%rot
  }

  
  ####PCA DO SPEC
  if (c("f4")%in%featype) {
    ncomps<-y$ncomps
    rot<-eigen(var(specdata))$vectors[,1:ncomps]
    loadspecdata<-rot
    pcaspecdata<-specdata%*%rot
  }
  
  ####WAVELET
  if (c("f7")%in%featype) {
  
    x<-as.matrix(wavCWT(1:LL,wavelet=y$wavelet,variance=y$variance))
    L0<-nrow(x)*ncol(x)
    cwtdata <- mat.or.vec(L0,nel)
  
    feacwt<-which(featype=="f7")
    wcwt<-c()
    for(i in feacwt){
      wcwt<-c(wcwt,y$W[[i]])
    }
    wcwt<-unique(floor((wcwt-1)/L0)+1)
  
    for (el in wcwt){
      cwtdata[,el]<-abs(as.vector(wavCWT(data[,el],wavelet=y$wavelet,variance=y$variance)))
    }
  }

  
  L<-y$L
  contwin<-0
  contW<-0
  contFea<-1
  features<-numeric(y$nfea)
  win<-y$win
  stat<-y$stat
  mintomax<-y$mintomax
  log<-y$log
  abs<-y$abs
  power<-y$power
  W<-y$W
  nel<-y$nch
  
  for (FEA in featype){
  
    if (FEA=="f1"){
      contW<-contW+1
      contwin<-contwin+1
      Ww<-W[[contW]]
      log[contwin]
      features[contFea:(contFea+length(Ww)-1)]<-.winPrac(data,win[contwin],stat[contwin],power[contwin],abs[contwin],log[contwin],mintomax[contwin],Ww,nel)
      contFea<-contFea+length(Ww)
      
    }
    if (FEA=="f2"){
      contW<-contW+1
      Ww<-W[[contW]]
      features[contFea:(contFea+length(Ww)-1)]<-abs(loaddata[Ww])
      contFea<-contFea+length(Ww)
      
    }
    if (FEA=="f3"){
      contW<-contW+1
      contwin<-contwin+1
      Ww<-W[[contW]]
      features[contFea:(contFea+length(Ww)-1)]<-.winPrac(specdata,win[contwin],
      stat[contwin],power[contwin],abs[contwin],log[contwin],mintomax[contwin],Ww,nel)
      contFea<-contFea+length(Ww)
    }
    if (FEA=="f4"){
      contW<-contW+1
      Ww<-W[[contW]]
      features[contFea:(contFea+length(Ww)-1)]<-abs(loadspecdata[Ww])
      contFea<-contFea+length(Ww)
    }
    if (FEA=="f5"){
      contW<-contW+1
      contwin<-contwin+1
      Ww<-W[[contW]]
      features[contFea:(contFea+length(Ww)-1)]<-.winPrac(pcadata,win[contwin],stat[contwin],power[contwin],abs[contwin],log[contwin],mintomax[contwin],Ww,ncomps)
      contFea<-contFea+length(Ww)
    
    }
    if (FEA=="f6"){
      contW<-contW+1
      contwin<-contwin+1
      Ww<-W[[contW]]
      dados<-.winPrac2(data,win[contwin],stat[contwin],power[contwin],abs[contwin],log[contwin],mintomax[contwin],nel,L)
      contwin<-contwin+1
      features[contFea:(contFea+length(Ww)-1)]<-.winPrac(dados,win[contwin],stat[contwin],power[contwin],abs[contwin],log[contwin],mintomax[contwin],Ww,nel)
      contFea<-contFea+length(Ww)
    }
    if (FEA=="f7"){
      contW<-contW+1
      Ww<-W[[contW]]
      features[contFea:(contFea+length(Ww)-1)]<-as.vector(cwtdata)[Ww]
      contFea<-contFea+length(Ww) 
    }
  } #for FEA in featype


  
  m<-y$model
  Ans<-predict(m, t(as.matrix(features)))
  prob <- predict(m, t(as.matrix(features)), probability = TRUE)
  prob<-attr(prob,"probabilities")

  NCL<-y$model$nclasses
  wc<-y$which.classes
  if(NCL==2){
    if (Ans=="A") pred<-c(prob[1],wc[1]) else pred<-c(prob[2],wc[2])
  }else{
    if (Ans=="A") pred<-c(prob[1],wc[1]) else if (Ans=="B") pred<-c(prob[2],wc[2]) else pred<-c(prob[3],wc[3])
  }
  
  return(pred)

}
