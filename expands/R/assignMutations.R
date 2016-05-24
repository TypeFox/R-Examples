assignMutations<-function( dm, finalSPs, max_PM=6, peakselection='localsum'){
  # P_CUTOFF=0.8;  #P_CUTOFF will cause high probabilities to be masked. these seem to be a symptom of a bad fit of the model. ##TODO; justify
  
  if (is.null(dim(finalSPs)) || nrow(finalSPs)==1) {
    spFreq = finalSPs[ "Mean Weighted"]
    precision=finalSPs["precision"]
  }  else {
    spFreq = finalSPs[, "Mean Weighted"]
    precision=finalSPs[1,"precision"]
  }
  spFreq=sort(spFreq);
  
  ##PM_B is the ploidy of the B-allele in SP; PM is the total ploidy in SP_cnv
  addCols=c("%maxP","SP","PM_B","SP_cnv","PM","PM_cnv","scenario");
  for (k in 1:length(addCols)){
    dm=.addColumn(dm,addCols[k],NA);
  }
  if (!any(colnames(dm)=="f")){
    dm=.addF(dm,  max_PM);
  }
  dm[,"SP"]=NA;  dm[,"SP_cnv"]=NA; ##delete any potentially existing SP info
  freq=c()
  for (sp in spFreq){
    freq=c(freq,seq(sp-precision/2,sp+precision/2,by=precision/20))
  }
  
  success=0;
  densities=matrix(matrix(NA,nrow(dm),length(freq)),nrow=nrow(dm),ncol=length(freq),dimnames=list(1:nrow(dm),freq));
  for(k in 1:nrow(dm)){
    ##Joined fit
    snvJ=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=max_PM),silent=TRUE)
    ##Separate fit
    f_CNV=NA; pm=NA;
    cnv=try(cellfrequency_pdf(NA,dm[k,"CN_Estimate"],NA,freq, max_PM=max_PM, snv_cnv_flag=2),silent=TRUE)
    snvSbeforeC=NULL;
    if(class(cnv)!="try-error" && any(!is.na(cnv$p)) && max(cnv$p, na.rm=T)>0){
      maxP_CNV = .selectCellFrequencyPeak(cnv$p,strategy=peakselection);#,max_p_threshold =P_CUTOFF)
      idx=which.min(abs(spFreq-freq[maxP_CNV$index]))
      #idx=which.min(abs(spFreq-freq[which.max(cnv$p)]))
      if (!isempty(idx)){
        f_CNV=spFreq[idx];
        # idx=which.min(abs(cnv$fit[,"f"]-f_CNV)); ##TODO: use deviation, not just cell-frequency to find index
        idx=which(abs(cnv$fit[,"f"]-f_CNV)<=precision/2); ##index of fits matching SP size
        idx=idx[which.min(cnv$fit[idx,"dev"])] ##index of fit of matching SP size with minimal residual (dev)
        if (!isempty(idx)){
          pm=cnv$fit[idx,"PM"]; 
          ##Fit under the assumption that SNV happened before CNV
          snvSbeforeC=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=NA, snv_cnv_flag=4, SP_cnv=f_CNV, PM_cnv=pm),silent=TRUE);
        }
      }
    }
    ##Max_PM is either 2 if SP with SNV is not a descendant of SP with CNV.... 
    snvS_noDesc=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq, max_PM=2, snv_cnv_flag=1),silent=TRUE)
    ##.. or Max_PM is PM of SP with CNV otherwise
    snvS_Desc=try(cellfrequency_pdf(dm[k,"AF_Tumor"],dm[k,"CN_Estimate"],dm[k,"PN_B"],freq[freq<=f_CNV+precision/2], max_PM=pm, snv_cnv_flag=1),silent=TRUE)
    ##Choose better solution between the two
    snvS=snvS_noDesc;
    if(class(snvS_Desc)!="try-error" && any(!is.na(snvS_Desc$p)) && max(snvS_Desc$p, na.rm=T)>0){
      if(class(snvS_noDesc)!="try-error" && any(!is.na(snvS_noDesc$p)) ){
        maxP_SnoDesc = .selectCellFrequencyPeak(snvS_noDesc$p,strategy=peakselection);#,max_p_threshold =P_CUTOFF)
        maxP_SDesc = .selectCellFrequencyPeak(snvS_Desc$p,strategy=peakselection);#,max_p_threshold =P_CUTOFF)
      }
      if(class(snvS_noDesc)=="try-error" || all(is.na(snvS_noDesc$p)) || maxP_SnoDesc$P<maxP_SDesc$P){ ##max(snvS_noDesc$p,na.rm=T)<max(snvS_Desc$p,na.rm=T)
        tmp=matrix(0,length(freq),1);        
        tmp[freq<=f_CNV+precision/2]=snvS_Desc$p;  snvS_Desc$p=tmp;  ##Complement to cover entire spFreq space
        snvS=snvS_Desc;
      }
    }
    
    maxP_J=list(index=-1, P=0); ##Maximum probability from joined fit
    maxP_S=list(index=-1, P=0); ##Maximum probability from separate fit
    maxP_SbeforeC=list(index=-1, P=0); ##Maximum probability from separate fit, under the assumption that CNV happened in descendant of SP with SNV
    
    if(class(snvJ)!="try-error" && any(!is.na(snvJ$p))){
      #maxP_J=max(snvJ$p,na.rm=T)
      maxP_J = .selectCellFrequencyPeak(snvJ$p,strategy=peakselection);#,max_p_threshold =P_CUTOFF)
    }
    if(class(snvS)!="try-error" && any(!is.na(snvS$p))){
      #maxP_S=max(snvS$p,na.rm=T)
      maxP_S = .selectCellFrequencyPeak(snvS$p,strategy=peakselection);#,max_p_threshold =P_CUTOFF)
    }
    if(!is.null(snvSbeforeC) && class(snvSbeforeC)!="try-error" && any(!is.na(snvSbeforeC$p))){
      #maxP_SbeforeC=max(snvSbeforeC$p,na.rm=T)
      maxP_SbeforeC = .selectCellFrequencyPeak(snvSbeforeC$p,strategy=peakselection);#,max_p_threshold =P_CUTOFF)
    }
    
    ##Skip if no solution found
    if (maxP_J$P==0 && maxP_S$P==0 && maxP_SbeforeC$P==0){
      dm[k,"SP"]=NA;
      next;
    }
    
    joinedFit=FALSE;
    if (class(snvJ)!="try-error" && (dm[k,"PN_B"]==1 || maxP_J$P>=max(maxP_S$P,maxP_SbeforeC$P,na.rm=T))){ ##SP carrying SNV and SP carrying CNV have same size, i.e. are identical:
      joinedFit=TRUE; ##LOH has to be associated with copy number variation --> SNV and CNV must be fit together
      snv=snvJ; 
      dm[k,"scenario"]=3;
    }else{ 
      if(maxP_S$P>=maxP_SbeforeC$P){
        snv=snvS;
        dm[k,"scenario"]=1;
      }else{
        snv=snvSbeforeC;
        dm[k,"scenario"]=4;
      }
    }
    
    ##Save end result:
    # idx=which.min(abs(spFreq-freq[which.max(snv$p)]))
    maxP_ = .selectCellFrequencyPeak(snv$p,strategy=peakselection);#,max_p_threshold =P_CUTOFF)
    idx=which.min(abs(spFreq-freq[maxP_$index]))
    dm[k,"SP"]=spFreq[idx];  
    idx=which(abs(snv$fit[,"f"]-dm[k,"SP"])<=precision/2); ##index of fits matching SP size
    idx=idx[which.min(snv$fit[idx,"dev"])] ##index of fit of matching SP size with minimal residual (dev)
    if(!isempty(idx)){
      dm[k,c("PM_B","PM")]=snv$fit[idx,c("PM_B","PM")]; ##(dm[k,"CN_Estimate"]*dm[k,"AF_Tumor"]-(1-dm[k,"SP"])*dm[k,"PN_B"])/dm[k,"SP"];  dm[k,"PM_B"]=max(1,dm[k,"PM_B"]);
    }
    if (!is.na(dm[k,"PM"]) && dm[k,"PM"]<0){
      dm[k,"PM"]=NA; ##PM can be -1 if obtained with snv_cnv_flag=1; TODO --> get NA directy for jar and remove this. 
    }
    #dm[k,"%maxP"]=max(snv$p,na.rm=T); 
    dm[k,"%maxP"]=maxP_$P
    densities[k,]=snv$p;
    
    if(joinedFit){
      dm[k,"SP_cnv"]=dm[k,"SP"];
      dm[k,"PM_cnv"]=dm[k,"PM"];
    }else if (!is.na(pm) && pm==2){
      dm[k,"SP_cnv"]=dm[k,"SP"];
      dm[k,"PM_cnv"]=pm; dm[k,"PM"]=pm;
      #dm[k,"PM"]=pm;
    }else{
      dm[k,"SP_cnv"]=f_CNV;
      dm[k,"PM_cnv"]=pm; 
      
      ##PM of SP does not have to be 2, because SP may also be a descendant of SP_cnv, i.e. the clone that acquired the CNV
      ##IF SP is larger than SP_cnv, then SP cannot be descandant of SP_cnv and therefor cannot harbor the CNV --> ploidy = 2
      ##IF SP is smaller than SP_cnv than SP may descend from SP_cnv: 
      ##-->Reject descendant hypothesis (ploidy of SP = 2) --> if 2 >= PM_B > PM 
      ##-->Accept descendant hypothesis (ploidy of SP = PM) --> if PM >= PM_B > 2 
      ##-->Irresolvable otherwise (ploidy of SP cannot be assigned) 
      if(!is.na(dm[k,"SP"]) && !is.na(dm[k,"SP_cnv"])){
        if(dm[k,"SP"]>dm[k,"SP_cnv"]){
          dm[k,"PM"]=2;
        }else{
          if(is.na(dm[k,"PM_B"]) || is.na(dm[k,"PM_cnv"])){
            dm[k,"PM"]=NA;
          }else if(dm[k,"PM_B"]>dm[k,"PM_cnv"] && dm[k,"PM_B"]<=2){
            dm[k,"PM"]=2;
          }else if(dm[k,"PM_B"]>2 && dm[k,"PM_B"]<=dm[k,"PM_cnv"]){
            dm[k,"PM"]=dm[k,"PM_cnv"];
          }else{
            dm[k,"PM"]=NA;
          }
        }
      }
    }
    
    success=success+1;
    if (mod(k,20)==0){
      print(paste("Processed", k, "out of ",nrow(dm),"SNVs --> success: ",
                  success,"/",k))
    }
  }
  
  dm[dm[,"%maxP"]==0,"SP"]=NA;
  
  ##Remove SPs to which no mutations were assigned
  toRm=c();
  for (j in 1:size(finalSPs,1)){
    if(is.null(dim(finalSPs))){
      idx=which(dm[,"SP"]==finalSPs["Mean Weighted"]);
      finalSPs["nMutations"]=length(idx);	  
    }else{
      idx=which(dm[,"SP"]==finalSPs[j,"Mean Weighted"]);
      finalSPs[j,"nMutations"]=length(idx);
      if(length(idx)==0){
        toRm=c(toRm,j);
      }
    }
  }
  if(length(toRm)>0){
    finalSPs=finalSPs[-1*toRm,, drop=FALSE];
  }
  output=list("dm"=dm,"finalSPs"=finalSPs);
  return(output);
}


.addF<-function (dm,  max_PM){
  #.jaddClassPath("ExPANydS.jar")
  .jinit(classpath="ExPANdS.jar")
  #javaImport(packages = "core.analysis.ngs.algorithms.*")
  dm=.addColumn(dm,"f",NA);
  snv_cnv_flag=3; ##co-occurrence assumption of SNV and CNV
  
  for (k in 1:nrow(dm)){
    expands <-try(.jnew("ExPANdS", as.double(dm[k,"AF_Tumor"]),as.double(dm[k,"CN_Estimate"]),
                        as.integer(dm[k,"PN_B"]),as.integer(max_PM)));
    if (class(expands)=="try-error"){
      print(expands);
      print(paste("At SNV ",k,": -->"));
      print(dm[k,]);
    }else{
      .jcall(expands,,"run",as.integer(snv_cnv_flag))
      dm[k,"f"]<-.jcall(expands,"D","getF");
    }
  }
  return(dm);
}

.selectCellFrequencyPeak<-function(probs, strategy='localsum',simple=TRUE){#,max_p_threshold = NULL){
  out=list(index=NA,P=NA)
  if(strategy=='maximum'){
    out=list(index=which.max(probs),P=max(probs,na.rm=T));
  }else if(strategy=='localsum'){
    # This function was devised to handle the observation that the probability distributions seem to favour LOH events even when the copy number is normal and generally seem to
    # favour a higher ploidy value. The peaks for the higher ploidy in the probability distributions are sharp but there is a clear local maximum often for the alternate ploidy
    # This function attempts to find the local maximum by calculating the sum of all probabilities in individual peaks and returning the index of the peak with the highest sum rather
    # than the maximum point probability, which is used in the default behaviour of Expands.
    
    #another feature available in this function is to assess the kurtosis of each peak. This is to handle another scenario that was observed in which some probability distributions contain
    #regions with very broad local maxima (blocky peaks). Such low kurtosis peaks seem to usually represent a poor quality fit and should ideally be removed. This function returns the kurtosis of the chosen peak
    #to help in determining if the fit is worth keeping. The feature is currently not being used and can sometimes lead to segfaults
    
    #"simple" mode skips all kurtosis calculations
    
    # if threshold is supplied, toss any peak region that contains a probability > this. Lower quality predictions seem to result from fits of the model with such high values. 
    
    #one known (or suspected) limitation of this function is that probabilities associated with different SPs are not considered separately and almost certainly should be. It's not clear how serious this problem may be
    
    #   if(!missing(max_p_threshold)){
    #     above_thresh = which(probs > max_p_threshold)
    #     probs[above_thresh] = 0
    #   }
    peak_regions = sign(diff(probs))
    
    ind = 1
    last_sign = 0
    last_ind = 0
    peakmax = c(0)
    peaksum = c(0)
    peakmax_ind = c(0)
    peaknum = 1
    peakvals=c()
    peak_kurtosis = c()
    for(i in peak_regions){
      if(i < 0){
        if(last_sign == -1){
          #same peak
          peaksum[peaknum] = peaksum[peaknum] + probs[ind]
          if(peakmax[peaknum]<probs[ind]){
            peakmax[peaknum] = probs[ind]
            
            peakmax_ind[peaknum] = ind
          }
          peakvals=c(peakvals,probs[ind])
          last_ind = ind
          ind = ind+1
        } else{
          #new peak
          peakmax[peaknum] = probs[ind]
          peaksum[peaknum] = peaksum[peaknum] + probs[ind]
          peakmax_ind[peaknum] = ind
          peakvals = c(probs[ind])
          last_ind = ind
          ind = ind+1
        }
        last_sign = i
      } else{
        if(last_sign < 0){
          if(!simple){
            kurtosis_last_peak = kurtosis(peakvals)
            peak_kurtosis[peaknum] = kurtosis_last_peak
          }
          peaknum=peaknum+1
          peaksum = c(peaksum,0)
        }
        last_sign = i
        last_ind = ind
        ind = ind+1
      }
    }
    if(!simple){
      kurtosis_last_peak = kurtosis(peakvals)
      peak_kurtosis[peaknum] = kurtosis_last_peak
    }
    best_peak_idx = which.max(peaksum)
    best_peak_max_idx = peakmax_ind[best_peak_idx]
    max_val_best_peak = peakmax[best_peak_idx]
    out=list(index=best_peak_max_idx,P=peaksum[best_peak_idx]);
    if(!simple){
      out$max_val_best_peak=max_val_best_peak
      out$kurtosis=peak_kurtosis[best_peak_idx]
    }
  }
  #   if (length(out)>5){
  #     x=unlist(out)
  #     out=list(index=as.numeric(x[grep('index',names(x))]),P=as.numeric(x[grep('P',names(x))]))
  #   }
  return(out);
}

