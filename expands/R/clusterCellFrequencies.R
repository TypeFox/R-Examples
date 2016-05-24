clusterCellFrequencies <- function(densities, precision, nrep=30, min_CellFreq=0.1){ ##, plotF=0
  #   if(plotF>0 && !require(rgl)){
  #   	plotF=0;
  #   	message("Plot supressed:  Package rgl required for 3D plot of subpopulation clusters. Load this package befor using this option.")
  #   }
  freq=as.numeric(colnames(densities));
  print(paste("Clustering ",nrow(densities),"probability distributions..."))
  cols=c("red","yellow","green","pink","magenta","cyan","lightblue","blue");
  
  count=0;##counting xtick for cluster plot
  clusterMethod="average";
  print(paste("Clustering agglomeration method:",clusterMethod))
  
  ##Cluster probabilities based on Kullback-Leibler divergence
  D = KLdiv(t(densities),eps=10^-18);
  idxx=which(apply(is.finite(D),1,all) );
  print(paste(nrow(densities)-length(idxx),"SNVs excluded due to non-finite pdfs") );
  D=D[idxx,idxx];densities=densities[idxx,];
  Z = hclust(as.dist(D),method = clusterMethod);
  TC = cutree(Z,k=round(sqrt(nrow(densities))));
  clIdx=unique(TC);
  print("Done");
  #figure("Name",paste("Cluster-Size_",R,sep=""));hist(T,length(clIdx))
  
  print("Filtering Clusters...");
  allSPs=list();tRep=nrep;
  while (nrep>0){
    if (mod(nrep,3)==0){
      print(paste(100*(tRep-nrep)/tRep, "% completed"))
    }
    ##test each cluster for significance --> printlay SPs
    spCols=c("Max Size","Mean Size","Mean Weighted","wilcTest","wilcTest_Mean","wilcTest_Sum","wilcTest_Kurtosis","kurtosis","kurtosisNoise","kurtosisMean","kurtosisNoiseMean","nMutations","precision","score");
    SPs <- matrix(nrow = length(clIdx), ncol = length(spCols), byrow = TRUE, dimnames = list(paste(c(1:length(clIdx))), spCols))
    
    for (k in 1:length(clIdx)){
      meanClPeak<-peak<-wMean<-score<-NA; ##init
      
      ia=which(TC==clIdx[k]);
      if (length(ia)<2){
        next;
      }
      clusterM=densities[ia,];
      ##Sort
      ia=apply(clusterM,1,which.max);
      ia=order(freq[ia]);
      clusterM=clusterM[ia,];
      
      ##Extend around maxima
      peak=.weightedMean(clusterM,freq);
      peakMin=peak-0.05;peakMax=peak + 0.05;
      idx=find(freq>=peakMin & freq<=peakMax);
      stdOk=which(apply(densities[,idx],1,std)>10^-5);
      if (length(stdOk)<5){
        next;
      }
      densitiesOk=densities[stdOk,];##remove densities with small std
      #obj=densityMclust(densitiesOk[,idx],G=2:2);
      #Tx=obj$classification;
      #Tx=Kmeans(densitiesOk[,idx],2,method="euclidean",nstart = 1,iter.max = 20)
      Tx=kmeans(densitiesOk[,idx], centers=2, nstart =10);
      Tx=Tx$cluster;
      if (length(unique(Tx))!=2 || length(which(Tx==1))<=1 ||
          length(which(Tx==2))<=1){
        next; #second cluster step unsuccessfull for this range;
      }
      
      cl1=densitiesOk[which(Tx==1),]; ia=which.max(apply(cl1,2,mean));
      cl2=densitiesOk[which(Tx==2),]; ib=which.max(apply(cl2,2,mean));
      meanCl=c(freq[ia],freq[ib]);
      peakIdx=which.min(abs(meanCl-peak));
      peakCl=densitiesOk[which(Tx==peakIdx),];
      meanClPeak=meanCl[peakIdx]; 
      
      wMean=.weightedMean(peakCl[,idx],freq[idx]);
      tryCatch({   
        ##find peak range of cluster
        maxCl=apply(peakCl,2,mean,na.rm=T);
        x=apply(peakCl[,idx],1,max,na.rm=T);  
        y=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,max,na.rm=T);  
        zz=wilcox.test(x,y,conf.level=0.99,alternative="greater");
        
        x2=apply(peakCl[,idx],1,mean,na.rm=T);  
        y2=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,mean,na.rm=T);  
        zz2=wilcox.test(x2,y2,conf.level=0.99,alternative="greater");
        
        x3=apply(peakCl[,idx],1,sum);
        y3=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,sum);              
        zz3=wilcox.test(x3,y3,conf.level=0.99,alternative="greater");
        kurt=apply(peakCl[,idx],1,kurtosis,na.rm=T);
        kurtNoise=apply(peakCl[,setdiff(c(1:length(freq)),idx)],1,kurtosis,na.rm=T);
        zzK=wilcox.test(kurt,kurtNoise,conf.level=0.99,alternative="greater");
        
        SPs[k,]=c(peak,meanClPeak,wMean, zz$p.value,zz2$p.value,zz3$p.value,zzK$p.value,max(kurt),max(kurtNoise),mean(kurt),mean(kurtNoise),nrow(peakCl),precision, score);
        ##calculate score
        score=SPs[k,"wilcTest"]+SPs[k,"wilcTest_Mean"]+SPs[k,"wilcTest_Sum"]+1/log(SPs[k,"nMutations"]);
        if(!is.na(SPs[k,"kurtosisNoiseMean"])){
          score=score+SPs[k,"kurtosisNoiseMean"]/500;
        }else{
          score=score+1;
        }
        SPs[k,"score"]=score;
        
        #         #plot option
        #         if (exists("plotF") && plotF>0 && nrep==tRep){
        #            col=cols[mod(k,length(cols))+1];
        #            count=.addTo3DPlot(count,clusterM,freq,col);
        #         }
      },error = function(e) {
        print(e);
      })
    }    
    ##collapse similar
    ia=order(SPs[,"Mean Weighted"],decreasing=T);
    SPs=SPs[ia,];
    SPs=.collapseSimilar(SPs,precision);
    ##print(paste("Found ",size(SPs,1),"SPs."));
    if (size(SPs,1)>0){
      allSPs[[length(allSPs)+1]]=SPs;
    }
    nrep=nrep-1;
  }
  
  if (length(allSPs)==0){
    return(NULL);
  }
  
  #   ##plot option
  #   if (plotF>0){
  #     title3d("",label);
  #   }
  
  robSPs=.chooseRobustSPs(allSPs,precision,min_CellFreq);
  SPs=.collapseSimilar(robSPs$SPs,precision);
  
  outcols=c("Mean Weighted","score","precision","nMutations");##printlay only these columns
  if (is.null(dim(SPs)) || nrow(SPs)==1){
    SPs=SPs[outcols];
  }else{
    SPs=SPs[,outcols];
  }
  
  print("Done.");
  return(SPs);
}

.weightedMean<-function(peakCl,freq){
  ##weighted mean
  wMean=0;sumWeight=sum(apply(peakCl,1,na.rm=T,max));
  for (pI in 1:nrow(peakCl)){
    maxIdx=which.max(peakCl[pI,]);
    wMean=wMean+(peakCl[pI,maxIdx]/sumWeight)*freq[maxIdx];
  }
  return(wMean);
}

.chooseRobustSPs <- function(allSPs,precision,min_CellFreq){
  ## input parameter SPs is a cell array with DataMatrix (DM) objects. Each row in
  ## DM contains the size and the p-value associated with a SP. DM is sorted
  ## in descending order of SP size.
  #count frequencies among predictions;
  freq=t(seq(min_CellFreq,1,by=precision));
  SPsizes=matrix(nrow = length(allSPs), ncol = length(freq),
                 dimnames = list(paste(c(1:length(allSPs))), freq))
  for (i in 1:length(allSPs)){
    SPs=allSPs[[i]];
    if(is.null(dim(SPs))){
      #       if(SPs["Mean Weighted"]>1){ ##Should no longer be necessary
      #         SPs["Mean Weighted"]=1; 
      #       }
      idx=which.min(abs(SPs["Mean Weighted"]-freq));
      SPsizes[i,idx]=SPs["score"];
    }else{
      #SPs[SPs[,"Mean Weighted"]>1,"Mean Weighted"]=1; ##Should no longer be necessary
      sps=SPs[,"Mean Weighted"];
      for (j in 1:length(sps)){
        idx=which.min(abs(sps[j]-freq));
        SPsizes[i,idx]=SPs[j,"score"];
      }
    }
  }
  keep=which(apply(!is.na(SPsizes),2,sum)>length(allSPs)*0.5);
  ##keep only best among recurrent SPs
  finalSPs=c();
  for (i in 1:length(keep)){
    sp_freq=freq[keep[i]];
    spI=which(!is.na(SPsizes[,keep[i]]));
    ia=which.min(SPsizes[spI,keep[i]]);
    SPs=allSPs[[as.numeric(spI[ia])]];
    if(is.null(dim(SPs))){
      SPs["Mean Weighted"]=sp_freq; ##standardize sp size
      finalSPs=rbind(finalSPs,SPs);
    }else{
      ia=which.min(abs(SPs[,"Mean Weighted"]-sp_freq));
      SPs[ia,"Mean Weighted"]=sp_freq; ##standardize sp size
      finalSPs=rbind(finalSPs,SPs[ia,]);
    }
  }
  output=list("SPs"=finalSPs,"spGrid"=SPsizes);
  return(output);
}


#.addTo3DPlot <- function(count,clusterM,freq,color,myPlot){
#X=(count+1):(count+nrow(clusterM))
#updatecount=count+nrow(clusterM);
#addV=TRUE;
#if(count==0){
#  addV=FALSE;
#}
#persp3d(as.numeric(X),as.numeric(freq),clusterM,col=color,aspect=c(1, 1, 0.5), add=addV,
#                  xlab="Mutation", ylab="cell-frequency", zlab="Probability");
#return(updatecount);
#}

.collapseSimilar <-function(SPs,precision){
  isNaNIdx=which(is.na(SPs[,"Mean Weighted"]));
  if (!isempty(isNaNIdx)){
    SPs=SPs[-isNaNIdx,];
  }
  if (size(SPs,1)<2){
    return(SPs);
  }
  spSize=unique(round(SPs[,"Mean Weighted"]*100)/100);
  for (n in 1:length(spSize)){
    idx=which(abs(SPs[,"Mean Weighted"]-spSize[n])<precision);
    if (length(idx)>1){
      ia=which.min(SPs[idx,"wilcTest"]);
      rmIdx=setdiff(idx,idx[ia]);
      SPs=SPs[-rmIdx,];
    }
    if (is.null(dim(SPs))){
      break;
    }
  }
  return(SPs);
}