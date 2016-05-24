.quantPriors<-function(x,quant){
  thresh<-quantile(x,quant)
  less<-x<=thresh
  great<-x>=thresh
  n<-length(x)
  prop<-c(sum(less)/n,sum(great)/n)
  mu<-c(mean(x[less]),mean(x[great]))
  sigma<-c(sd(x[less]),sd(x[great]))
  return(list(prop,mu,sigma))
}
.tooClose<-function(mod){
  parms<-do.call("rbind",list(mod$mu,mod$sigma,mod$lambda))
  big.var<-which.max(parms[2,])
  min.var<-which.min(parms[2,])
  uniCheck<-abs((parms[1,big.var]-parms[1,min.var])) < parms[2,big.var]
  tooSmall<-any(parms[3,]<.01)
  return(uniCheck|tooSmall)}
.noiseSeg<-function(rtsnr,step=.05){
  print("Running Expectation Maximisation on Robust TSNR for noise mask creation")
  resamp<-quantile(rtsnr,seq(0,.98,length.out = 1000))
  
  props<-seq(.01,.95,step)
  vals<-list()
  for(i in 1:length(props)){
    vals[[i]]<-.quantPriors(resamp,props[i])
  }
  
  mods<-list()
  for(i in 1:length(props)){mods[[i]]<-gmm(x = resamp,k = 2,imeans = vals[[i]][[2]],isd =vals[[i]][[3]],ilambda = vals[[i]][[1]])}
  
  tc<-sapply(mods,.tooClose)
  if(sum(tc)==length(vals)){return("a Mixture model may not be appropriate. Therefore the threshold is defined using a single gaussian")}
  possibleMods<-mods[tc==FALSE]
  ll<-sapply(possibleMods,function(x){x$LL})
  mus<-sapply(possibleMods,function(x){sort(x$mu)})
  muDiff<-abs(mus[1,]-mus[2,])
  modInd<-which.min(ll*muDiff)
  return(possibleMods[[modInd]])
}
.nnInterp<-function(mask.arr,n,method=1){
  #get dimensions
  dim.arr<-dim(mask.arr)
  dm<-dilate(mask.arr,k = 1)
  #m<-dm>0
  inds<-t(which(dm>0,arr.ind = TRUE))
  edge<-which((inds[1,]==dim.arr[1]|inds[1,]==1) | (inds[2,]==dim.arr[2]|inds[2,]==1) | (inds[3,]==dim.arr[3]|inds[3,]==1) ,arr.ind = TRUE)
  inds.edge<-inds[,edge]
  if(ncol(inds.edge)>0){inds<-inds[,-edge]}
  rownames(inds)<-c("dim1","dim2","dim3")
  neighbours<-t(expand.grid(c(n:-n),c(n:-n),c(n:-n)))
  if(method==2){neighbours<-neighbours[,neighbours[3,]==0]
  } else if(method==3){neighbours<-neighbours[,neighbours[2,]==0]
  } else if(method==4){neighbours<-neighbours[,neighbours[1,]==0]}
  
  inds.neighbours<-vector("list",ncol(neighbours))
  for(i in 1:length(inds.neighbours)){
    inds.neighbours[[i]]<-t(inds+neighbours[,i])
  }
  inds.edge.neighbours<-vector("list",ncol(neighbours))
  for(i in 1:length(inds.edge.neighbours)){
    inds.edge.neighbours[[i]]<-t(inds.edge+neighbours[,i])
  }
  ##bounds check
  if(ncol(inds.edge)>0){
    for(i in 1 :length(inds.edge.neighbours)){inds.edge.neighbours[[i]][inds.edge.neighbours[[i]][,1]>dim.arr[1],1]<-dim.arr[1]}
    for(i in 1 :length(inds.edge.neighbours)){inds.edge.neighbours[[i]][inds.edge.neighbours[[i]][,2]>dim.arr[2],2]<-dim.arr[2]}
    for(i in 1 :length(inds.edge.neighbours)){inds.edge.neighbours[[i]][inds.edge.neighbours[[i]][,3]>dim.arr[3],3]<-dim.arr[3]}
    for(i in 1 :length(inds.edge.neighbours)){inds.edge.neighbours[[i]][inds.edge.neighbours[[i]]<1]<-1}
  }
  neighbour.val<-vector("list",ncol(neighbours))
  for(i in 1:length(neighbour.val)){
    if(ncol(inds.edge)>0){  
      neighbour.val[[i]]<-c(mask.arr[inds.neighbours[[i]][,,drop=FALSE]],mask.arr[inds.edge.neighbours[[i]][,,drop=FALSE]])
    }else{neighbour.val[[i]]<-mask.arr[inds.neighbours[[i]][,,drop=FALSE]]}
  }
  neighbour.mat<-do.call("cbind",neighbour.val)
  sums.neighbours<-rowSums(neighbour.mat)
  interp.mask.vec<-matrix(ifelse(sums.neighbours>floor(ncol(neighbours)/2),1,0),nrow=1)
  if(ncol(inds.edge)>0){
    mask.arr[t(cbind(inds,inds.edge))]<-interp.mask.vec
  }else{mask.arr[t(inds)]<-interp.mask.vec}
  return(mask.arr)
}

fiach <-function(input,t,tr,rp=NULL,maxgap=1,freq=128,nMads = 1.96){
  ################################################################
  ########### DIA DHUIT A CHARA, CHONAS ATA TU? ##################
  ########### BHUEL NA BI AG CAOINE.            ##################
  ###########  !!!!!TA FIACH ANSEO ANOIS!!!!!   ##################
  ################################################################
  
  ######################################
  ######## ARGUMENT CHECKS #############
  ######################################
  if(nMads<0){stop("nMads must be a postive real number")}
  if(t<0){stop("t must be 0 or a postive real number")}
  if(tr<0){stop("tr must be a postive real number")}
  if(maxgap<1){stop("maxgap must be greater than 1 ")}
  if(freq<0){stop("freq must be a postive real number")}
  if(!is.null(rp)){
  if(!is.character(rp)){stop("rp should be a character strings")}
  if(length(rp)>1){stop("Only one rp file should be specified")}
  if(!file.exists(rp)){stop("The specified rp file does not exist")}
  }
  if(length(input)<1){stop("At least one functional file must be specified")}
  if(!is.character(input)){stop("input should be character strings")}
  exists<-file.exists(input)
  if(!all(exists)){stop("At least one of the specified functional files does not exist")}
  #####################################
  ######## DATA READ ##################
  #####################################
  if(length(input)==1){four.d<-TRUE}else{four.d<-FALSE}
  data<-readNii(input)
  temp<-RNiftyReg::updateNifti(data[,,,1],template = data)
  outType<-getDatatype(input = input,type = "RNiftyReg")
  print("Data is Read") 
  #####################################
  ######### MATRIX CREATION  ##########
  #####################################                                                                               
  time<-dim(data)[4]                                                                            
  mat<-zeroNa(arrMat(data))
  nvox<-cumprod(dim(data)[1:3])[3]
  #####################################
  ###### FILE MANAGEMENT ##############
  #####################################
  dir<-dirname(input)[1]
  name<-basename(input)
  gzipped<-grepl("[.]gz$",name[1])
  if(gzipped){
    file.base<-substr(name,1,(nchar(name)-7))
  }else{file.base<-substr(name,1,(nchar(name)-4))}
  output.file.name<-paste(dir,"/filt_",file.base,sep="")                                 
  if(gzipped){output.file.name<-paste(dir,"/filt_",file.base,".nii.gz",sep="")}
  fold<-paste(dir,"_fiach_diagnostics",sep="")            
  suppressWarnings(dir.create(fold))               
  print("Output Folder Created")
  #####################################
  ######## MASK CREATION ##############
  #####################################
  meds<-colMedian(mat)                                                                                             
  mask.mat<-kmeansMask(meds)
  vf<-sum(mask.mat)/length(mask.mat)
  bss0<-(mean(meds[mask.mat==0])-mean(meds))^2*(sum(mask.mat==0))
  bss1<-(mean(meds[mask.mat==1])-mean(meds))^2*(sum(mask.mat==1))
  bss<-bss0+bss1
  tss<-sum((meds-mean(meds))^2)
  fit<-bss/tss
  if(vf<.25 & fit< .8 ){
    print("As the k-means clustering did not have a good fit the mask was construced using quantiles. Is there a large receive field bias in your data?")
    mask.mat<-quantMask(meds)}
  mask.arr<-matArr(mask.mat,dim=c(dim(data)[1:3],1)) 
  small.brain<-mat[,mask.mat==1]
  print("Brain Extraction Completed")
  #####################################
  ####### HIGH PASS FILTERING #########
  #####################################
  small.hp.mat<-highPass(small.brain,freq=freq,tr=tr) ## CHANGE for PSEUDO + HighBasis                                                        
  hp.mat<-mat
  hp.mat[,mask.mat==1]<-small.hp.mat
  hp.mat[,mask.mat==0]<-0
  hp.arr<-matArr(hp.mat,dim=dim(data))                                                       
  print(paste("High pass filtered at",freq, "seconds"))
  #####################################
  ########## TSNR MAPS ################
  ####################################
  brain.meds<-meds[mask.mat==1]                                                              
  brain.mads<-colMad(small.hp.mat)                                                               
  brain.robust.tsnr<-brain.meds/brain.mads                                                      
  output.brain.robust.tsnr<-numeric(nvox)                                                         
  output.brain.robust.tsnr[which(mask.mat==FALSE)]<-0                                           
  output.brain.robust.tsnr[which(mask.mat==TRUE)]<-brain.robust.tsnr                                    
  robust.brain.tsnr.arr<-matArr(output.brain.robust.tsnr,dim=c(dim(data)[1:3],1))             
  output.brain.meds<-numeric(nvox)                                                              
  output.brain.meds[which(mask.mat==FALSE)]<-0                                                  
  output.brain.meds[which(mask.mat==TRUE)]<-brain.meds                                          
  brain.meds.arr<-matArr(output.brain.meds,dim=c(dim(data)[1:3],1)) 
  print("Robust TSNR calcualted")                   
  #################################################
  ########## MIXTURE MODELLING ####################
  #################################################
  print("Estimating Mixture Model")                        
  mixmdl<-.noiseSeg(brain.robust.tsnr)                     
  if(is.character(mixmdl)){
  print(mixmdl)
   thresh<-quantile(brain.robust.tsnr,.05)
   segment.file<-paste(fold,"/segmentation.pdf",sep="")
   pdf(segment.file,width=11.69,height=8.27)
   hist(brain.robust.tsnr)
   abline(v=thresh)
   dev.off()}else{
   maxProp<-which.max(mixmdl$lambda)
   minProp<-which.min(mixmdl$lambda)
   thresh<-qnorm(.05,mixmdl$mu[maxProp],mixmdl$sigma[maxProp])
  
   segment.file<-paste(fold,"/segmentation.pdf",sep="")
   pdf(segment.file,width=11.69,height=8.27)
   plot(mixmdl)
   abline(v=thresh)
   dev.off()
   }
  print("Mixture Model is Estimated")
  ###################################################
  ####### CREATE NOISE MASK AND MODEL IT ############
  ###################################################
  print("Creating Noise Mask")
  noise.mask<-brain.robust.tsnr
  noise.mask[brain.robust.tsnr<thresh]<-1
  noise.mask[brain.robust.tsnr>thresh]<-0
  output.brain.noise<-integer(nvox)                                                                  
  output.brain.noise[which(mask.mat==0)]<-0                                              
  output.brain.noise[which(mask.mat==1)]<-noise.mask                                     
  brain.noise.arr<-matArr(output.brain.noise,dim=c(dim(data)[1:3]))
  
  cubenn.brain.noise.arr<-.nnInterp(brain.noise.arr,n=1,method=1)  
  in1nn.brain.noise.arr<-.nnInterp(brain.noise.arr,n=1,method=3)  
  in2nn.brain.noise.arr<-.nnInterp(brain.noise.arr,n=1,method=4)
  inplane.brain.noise.arr<-.nnInterp(brain.noise.arr,n=1,method=2)
  total<-cubenn.brain.noise.arr+in1nn.brain.noise.arr+in2nn.brain.noise.arr+inplane.brain.noise.arr
  total<-ifelse(total>=1,1,0)
  
  noise.mask.vec<-as.vector(total)
  small.noise.mask<-noise.mask.vec[mask.mat==1]

  noise<-small.hp.mat[,small.noise.mask==1]
 
  scale.noise<-scale(noise)
  non.zero.scale.noise<-scale.noise[,colSums(noise)>0]
  sds.noise<-apply(noise,2,sd)
  svd.thresh<-sds.noise<=quantile(sds.noise,.98)
  
  try(noise.mod<-prcomp(non.zero.scale.noise), silent=TRUE)
  if(exists("noise.mod")==FALSE){
    noise.mod<-prcomp(non.zero.scale.noise[,svd.thresh])
  }
    noise.regs6<-noise.mod$x[,1:6]
  print("Noise Mask Created and Regressors Computed")
  #####################################
  ########## REPLACE WITH NA ##########
  #####################################
  print("Identifying Bad Data")
  na.brain<-badData(small.hp.mat, meds=brain.meds, mads=brain.mads, nMads=nMads, t=t)
  print("Found Bad Data")
  #####################################
  ###### SCRUB FIRST AND LAST #########
  #####################################
  print("Correcting Bad Data")
  na.1<-which(is.na(na.brain[1,])==TRUE)
  na.brain[1,na.1]<-brain.meds[na.1]
  na.time<-which(is.na(na.brain[time,])==TRUE)
  na.brain[time,na.time]<-brain.meds[na.time]
  #####################################
  ######### SPLINE the NAs ############
  #####################################
  filt.brain<-naSpline(na.brain,maxgap=maxgap)                      
  #####################################
  ####### Identify ANNOYING VOXELS ####
  #####################################
  scrub.inds<-which(is.na(filt.brain)==TRUE,arr.ind=TRUE)
  #####################################
  ########## SCRUB ANNOYING ###########
  #####################################
  scrubbed.brain<-filt.brain
  scrubbed.brain[scrub.inds]<-brain.meds[scrub.inds[,2]]
  output.brain.mat<-hp.mat                                                                    
  output.brain.mat[,which(mask.mat==FALSE)]<-0                                           
  output.brain.mat[,which(mask.mat==TRUE)]<-scrubbed.brain                               
  mov.arr<-matArr(output.brain.mat, dim(data))  
  print("Data Corrected... File Writing Begins")
  ####################################
  ######## 3d Writing ################
  ####################################
  if (four.d==FALSE){
    for(i in 1:length(input)){
      out<-RNiftyReg::updateNifti(mov.arr[,,,i],template=temp)
       RNiftyReg::writeNifti(out, file=output.file.name[i],datatype = outType)
    }
  }
  #####################################
  ######## DATA QUALITY METRIC ########
  #####################################
  metric<-matrix(ncol=4,nrow=1)
  metric[,1]<-sum(is.na(na.brain))
  metric[,2]<-(time*(ncol(small.brain)-sum(noise.mask.vec)))
  metric[,3]<-metric[,1]/metric[,2]*100
  metric[,4]<-mixmdl$mu[maxProp]
  metric<-metric[,3:4,drop=FALSE]
  colnames(metric)<-c("% Data Changed","Peak TSNR")
  write.table(metric,file=paste(fold,"/metrics.txt",sep=""),col.names=TRUE,row.names=FALSE)
  #####################################
  ######### RP MANIPULATION ###########
  #####################################
  if(is.null(rp)==FALSE){
  rp<-read.table(rp)
  if(nrow(rp)!=time){stop("The number of functional volumes does not equal the number of timepoints in the supplied rp file")}
  fd<-fd(rp)
  noise.basis6<-cbind(rp,noise.regs6)
  fdNoise<-cbind(fd,noise.regs6)
  scaleHP<-scale(scrubbed.brain,center = brain.mads,scale = brain.mads)
  gs<-rowMeans(scaleHP)
  write.table(gs,file=paste(dir,"/gs.txt",sep=""),col.names=FALSE,row.names=FALSE)
  write.table(fdNoise,file=paste(dir,"/fd_noise.txt",sep=""),col.names=FALSE,row.names=FALSE)
  write.table(noise.basis6,file=paste(dir,"/noise_basis6.txt",sep=""),col.names=FALSE,row.names=FALSE)
  }else{
    noise.basis6<-noise.regs6
    write.table(noise.basis6,file=paste(dir,"/noise_basis6.txt",sep=""),col.names=FALSE,row.names=FALSE)
    }
  
  #####################################
  ######## FILE WRITING ###############
  #####################################
  threeDim<-c(3,dim(data)[1:3],rep(1,4))
  if(four.d==TRUE){
    out<-RNiftyReg::updateNifti(mov.arr,template = temp)
    RNiftyReg::writeNifti(out, file=output.file.name,datatype = outType)
    } 
  
  RNiftyReg::writeNifti(mask.arr,paste(fold,"/mask",sep=""),template=temp,datatype ="char") 
  RNiftyReg::writeNifti(robust.brain.tsnr.arr,paste(fold,"/rtsnr",sep=""), template=temp,datatype = "float")
  RNiftyReg::writeNifti(brain.meds.arr,paste(fold,"/median",sep=""), template=temp,datatype = outType) 
  #####################################
  #### MEMORY MANAGEMENT ##############
  #####################################
  print("Memory Cleanup") 
  rm(list=ls())
  garbage<-gc()
  print("End")
}
