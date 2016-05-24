########################
############ 2. boyce index
########################    

#### functions calculating Boyce index (Hirzel et al. 2006) by Blaise Petitpierre & Frank Breiner(28.06.2013)

#### internal function calculating predicted-to-expected ratio for each class-interval

boycei<-function(interval,obs,fit){
  
  fit.bin<-fit
  obs.bin<-obs
  fit.bin[fit[]>=interval[1]&fit[]<=interval[2]]<-"i";fit.bin[fit.bin!="i"]<-0
  obs.bin[obs[]>=interval[1]&obs[]<=interval[2]]<-"i";obs.bin[obs.bin!="i"]<-0
  
  pi<-length(which(obs.bin=="i"))/length(obs)
  ei<-length(which(fit.bin=="i"))/length(fit.bin)
  fi<-pi/ei
  
  return(fi)
}

#### Calculating Boyce index as in Hirzel et al. 2006
#### fit: vector containing the suitability values of the study area (or extent used for the validation)
#### obs: vector containing the suitability values of the validation points
#### nclass : number of classes or vector with classes threshold. 
####          If nclass=0, Boyce index is calculated with a moving window (see next parameters)
#### windows.w : width of the moving window (by default 1/10 of the suitability range)
#### res : resolution of the moving window (by default 100 focals)
#### PEplot : if True, plot the predicted to expected ratio along the suitability class

ecospat.boyce<-function(fit,obs,nclass=0,window.w="default",res=100,PEplot=T){
  if(window.w=="default"){window.w<-(max(fit)-min(fit))/10}
  interval<-c(min(fit),max(fit))
  mini<-interval[1];maxi<-interval[2]
  
  if(nclass==0){
    vec.mov<- seq(from = mini, to = maxi - window.w, by = (maxi-mini- window.w)/res)
    vec.mov[res+1]<-vec.mov[res+1]+1
    interval<-cbind(vec.mov,vec.mov+window.w)
  }else if (length(nclass)>1){
    vec.mov<-c(mini,nclass)
    interval<-cbind(vec.mov,c(vec.mov[-1],maxi))
  }else if(nclass>0 & length(nclass)<2){
    vec.mov<-seq(from = mini, to = maxi, by = (maxi-mini)/nclass)
  }
  
  
  f<-apply(interval,1,boycei,obs,fit)
  
  if(length(f[which(f!="NaN")])<=2){b<-NA
    }else{
    b<-cor(f[which(f!="NaN")],vec.mov[which(f!="NaN")],method="spearman")}
  
  ID <- seq(1:(length(vec.mov)))
  HS <- apply(interval,1,sum)/2
  
  if (PEplot==T)plot((apply(interval[which(f!="NaN"),],1,sum)/2),f[which(f!="NaN")],
                     xlab="Habitat suitability",ylab="Predicted/Expected ratio")
  
  results<-list(F.ratio=f, Spearman.cor=round(b,3), HS=HS, ID=ID)
  
  return(results)
}
