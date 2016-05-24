BPEC.MCMC <- function(RawSeqs,CoordsLocs,MaxMig,iter,ds,PostSamples=0,dims=-1)
{
  if(dims==-1){
    dims=round(length(which((round(CoordsLocs,digits=0)==CoordsLocs)==FALSE))/nrow(CoordsLocs),digits=0)
  }
  if(dims!=-1)
  {
    dimscal=round(length(which((round(CoordsLocs,digits=0)==CoordsLocs)==FALSE))/nrow(CoordsLocs),digits=0)
    if(dims!=dimscal)
    {
      writeLines("\nThe dataset appears to have some integer-valued environmental/phenotypic entries.\nIf this is not correct, check the value you have given for dims.")
    }
  }
  if(PostSamples>iter/10)
  {
    writeLines("You cannot have PostSamples greater than iter/10, so it will be changed to iter/10.")
    PostSamples=iter/10
  }
  if(iter%%10!=0)
  {
    writeLines("iter needs to be a multiple of 10. The MCMC will not be run. ")
    return(-1)
  }
  Coordinates = CoordsLocs[,1:dims]
  CoordsLocs[is.na(CoordsLocs)]=-1
  if(ncol(CoordsLocs)==dims)
  {
    CoordsLocs=cbind(CoordsLocs,seq(1,nrow(CoordsLocs)))
    writeLines("The dataset did not contain a list of haplotypes per row,\nso the program will assume that row 1 corresponds to haplotype 1, row 2 to haplotype 2 etc. ")
  }
  EndRow=1:nrow(CoordsLocs)
  EndRow[]=-2
  Coords = CoordsLocs[,1:2]
  CoordsLocs=cbind(CoordsLocs,EndRow)
  locno=as.numeric(nrow(CoordsLocs))
  MaxLoc = as.numeric(ncol(CoordsLocs))
  CoordsLocs=as.numeric(as.vector(t(CoordsLocs)))
  
  SeqNames=as.numeric(gsub("[^0-9]","",names(RawSeqs))  )
  SeqCount=as.numeric(length(RawSeqs))  
  SeqLength=as.numeric(length(RawSeqs[[1]]))  
  #  UnSeqs=RawSeqs
  #  UnCount=SeqCount
  
  
  z= unlist(RawSeqs)
  z= matrix(z, nrow=length(RawSeqs[[1]]), ncol=length(RawSeqs))
  z= t(z)
  
  DelCol=NULL
  for(i in 1:ncol(z))
  {
    if(length(unique(z[,i]))==1)
    {    
      DelCol=c(DelCol,i) 
    }  
  }
  rm(z)
  if(is.null(DelCol)==FALSE)
  {
    for(i in 1:length(RawSeqs))
    {
      RawSeqs[[i]]=RawSeqs[[i]][-DelCol]
    }
  }
  
  # for(i in 1:504){
  # 
  # if(length(unique(z[,i]))>1)
  # {
  # print(i)
  # }
  # }
  SeqLength=as.numeric(length(RawSeqs[[1]]))  
  Seqs=array(NA,dim=(SeqCount*SeqLength))
  counter=1
  for (i in 1:SeqCount)
  {
    for(j in 1:SeqLength){
      Seqs[counter]=RawSeqs[[i]][j]
      counter=counter+1
    }
  }
  
  # rm(RawSeqs)
  
  count=SeqCount
  seeds=2
  ancestral=1
  
  MCMCout = .C("BPEC",ModeInitial=1,SeqR=Seqs,CoordsLocsR=CoordsLocs,CoordsDimsR=as.numeric(dims),SeqsFileR=SeqNames,SeqCountR=count,SeqLengthR=SeqLength,LocNoR=locno,MaxLocR=MaxLoc,maxmigR=MaxMig,seedsR=seeds,iterR=iter,dsR=ds,ancestralR=ancestral,SampleMeansR=numeric(2),SampleCovsR=numeric(2),SampleIndicesR=numeric(2),SampleClusterCodaR=numeric(2),SampleRootCodaR=numeric(2),PostSamplesR=PostSamples,levelsR=numeric(2),cladoR=numeric(2),EdgeTotalProbR=numeric(2),NoSamplesR=numeric(10*count+100),ClusterProbsR=numeric(2),countR=numeric(2),MigPMigProbsR=numeric(MaxMig+1),RootProbsR=numeric(2),RootLocProbsR=numeric(locno),MCMCparamsR=numeric(8),SeqLabelsR=numeric(10*count+100),NSeqR=numeric(10),errorcodeR=numeric(1),PACKAGE="BPEC")
  
  MCMCout$countR=max(c(count,MCMCout$countR[1]))
  
  if(MCMCout$errorcodeR[1]!=1)
  {
    MCMCout = .C("BPEC",ModeInitial=0,SeqR=Seqs,CoordsLocsR=CoordsLocs,CoordsDimsR=as.numeric(dims),SeqsFileR=SeqNames,SeqCountR=SeqCount,SeqLengthR=SeqLength,LocNoR=locno,MaxLocR=MaxLoc,maxmigR=MaxMig,seedsR=seeds,iterR=iter,dsR=ds,ancestralR=ancestral,SampleMeansR=numeric((PostSamples+1)*dims*(MaxMig+1)*seeds),SampleCovsR=numeric((PostSamples)*dims*dims*(MaxMig+1)*seeds),SampleIndicesR=numeric(sum((MCMCout$NoSamplesR+1))*2*(PostSamples)),SampleClusterCodaR=numeric((PostSamples)*(dims+4+(dims-2))*(MaxMig+1)*seeds),SampleRootCodaR=numeric(PostSamples*seeds),PostSamplesR=PostSamples,levelsR=numeric(MCMCout$countR+1),cladoR=numeric((MCMCout$countR+1)*(MCMCout$countR+1)),EdgeTotalProbR=numeric((MCMCout$countR+1)*(MCMCout$countR+1)),NoSamplesR=MCMCout$NoSamplesR,ClusterProbsR=numeric((MCMCout$countR+1)*(MaxMig+1)),countR=MCMCout$countR,MigProbsR=numeric(MaxMig+1),RootProbsR=numeric(2*MCMCout$countR),RootLocProbsR=numeric(locno),MCMCparamsR=numeric(8),SeqLabelsR=numeric(MCMCout$countR+1),NSeqR=as.numeric(MCMCout$countR),errorcodeR=numeric(1),PACKAGE="BPEC")
    
    names(MCMCout$MCMCparamsR)=c('PsiPrior','CentralSieve','TrialAdd','AveClusterWeight','ClustAccRate','RootAccRate','MeanGamma','MeanW')
    
    MCMCout$CoordsLocsR = t(array(MCMCout$CoordsLocsR,dim=c(length(MCMCout$CoordsLocsR)/MCMCout$LocNoR,MCMCout$LocNoR)))
    MCMCout$CoordsLocsR = MCMCout$CoordsLocsR[,-ncol(MCMCout$CoordsLocsR)]
    
    MCMCout$ClusterProbsR = t(array(MCMCout$ClusterProbsR,dim=c(MaxMig+1,MCMCout$countR+1)))
    MCMCout$countR = MCMCout$countR[1]
    
    #  MCMCout$rootfreqsR = MCMCout$rootfreqsR[1:MCMCout$countR]
    #MCMCout$clustersR=MCMCout$clustersR[1:sum(MCMCout$NoSamplesR)]
    MCMCout$levelsR = MCMCout$levelsR[1:MCMCout$countR]
    MCMCout$ClusterProbsR = MCMCout$ClusterProbsR[1:MCMCout$countR,]
    MCMCout$cladoR = MCMCout$cladoR[1:(MCMCout$countR*MCMCout$countR)]
    MCMCout$EdgeTotalProbR = MCMCout$EdgeTotalProbR[1:(MCMCout$countR*MCMCout$countR)]
    MCMCout$EdgeTotalProbR =array(MCMCout$EdgeTotalProbR,dim=c(MCMCout$countR,MCMCout$countR))
    
    #  MCMCout$RootProbsR = MCMCout$RootProbsR[1:MCMCout$countR]
    
    if(length(MCMCout$SeqsFile)<max(c(MCMCout$SeqLabelsR,MCMCout$SeqLabelsR)))
    {
      MCMCout$SeqsFile=c(MCMCout$SeqsFile,rep(0,max(c(MCMCout$SeqLabelsR,MCMCout$SeqLabelsR))-length(MCMCout$SeqsFile)))
    }
    
    MCMCout$SampleMeansR=is.finite(MCMCout$SampleMeansR)*MCMCout$SampleMeansR
    MCMCout$SampleCovsR=is.finite(MCMCout$SampleCovsR)*MCMCout$SampleCovsR
    
    MCMCout$SampleMeansR=array(MCMCout$SampleMeansR,dim=c(dims,MaxMig+1,seeds*PostSamples))
    MCMCout$SampleClusterCodaR=array(MCMCout$SampleClusterCodaR,dim=c(dims+2 + dims,MaxMig+1,seeds*PostSamples))
    MCMCout$SampleIndicesR=array(MCMCout$SampleIndicesR,dim=c(sum(MCMCout$NoSamplesR),seeds*PostSamples))
    
    MCMCout$SampleCovsR=array(MCMCout$SampleCovsR,dim=c(dims,dims,MaxMig+1,seeds*PostSamples))
    
    
    #finally, re-order the posterior samples using Papastamoulis (2010) label-switching
    
    # ActualMaxMig = max(MCMCout$SampleIndicesR)
    # mcmc.pars = array(0,dim=c(2*PostSamples,ActualMaxMig,dims+dims*dims))
    # z = array(0,dim=c(2*PostSamples,sum(MCMCout$NoSamplesR)))
    # 
    # for (i in 1:(2*PostSamples))
    # {
    #   for (j in 1:dims)
    #   {
    #     mcmc.pars[i,1:ActualMaxMig,j] <- MCMCout$SampleMeansR[j,1:ActualMaxMig,i] 
    #   }
    #   
    #   for (j in 1:dims)
    #   {
    #     for(l in 1:dims)
    #     {
    #       mcmc.pars[i,1:ActualMaxMig,dims+(j-1)*dims+l] <- MCMCout$SampleCovsR[j,l,1:ActualMaxMig,i]   
    #     }
    #   }
    #   z[i,] <- MCMCout$SampleIndicesR[,i]
    # }
    #
    # writeLines("\nProcessing posterior samples...")
    #
    # run<-ecr.iterative.1(z, ActualMaxMig)
    
    #reorder the MCMC output according to this method:
    # reordered.mcmc<-permute.mcmc(mcmc.pars,run$permutations)
    
    # newz = MCMCout$SampleIndicesR
    # for (i in 1:2*PostSamples)
    #{
    #  for (j in 1:dims)
    #  {
    #      MCMCout$SampleMeansR[j,1:ActualMaxMig,i] = reordered.mcmc$output[i,1:ActualMaxMig,j] 
    #  }
    
    #   for (j in 1:dims)
    #   {
    #     for(l in 1:dims)
    #     {
    #         MCMCout$SampleCovsR[j,l,1:ActualMaxMig,i]   = reordered.mcmc$output[i,1:ActualMaxMig,dims+(j-1)*dims+l] 
    #     }
    #   }
    
    #  for(j in 1:max(z))
    #  {      
    #   newz[MCMCout$SampleIndicesR[,i]==j,i] <- which(run$permutations[i,]==j)
    # }
    #}
    #MCMCout$SampleIndicesR=newz
    
    
    MaxVar = 0
    ChainMeans = array(0, dim=c(dim(MCMCout$SampleMeansR)[1],dim(MCMCout$SampleMeansR)[2],2))
    flag = 0
    for(i in 1:(MaxMig+1))
    {
      for(j in 1:dims)
      {
        for(l in 1:2)
        {
          ChainMeans[j,i,l]=mean(MCMCout$SampleMeansR[j,i,(1+(l-1)*PostSamples):(l*PostSamples)],na.rm=TRUE)
        }
        if(sum(is.na(ChainMeans[j,i,]))==0)
        {
          if(abs(ChainMeans[j,i,1]-ChainMeans[j,i,2])>0.05*(max(Coordinates[,j])-min(Coordinates[,j])))
          {
            writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations")
            flag=1
            break
          }
        }
      }
      if(flag==1)
      {
        break
      }        
    }
    MCMCout$MCChainMeansR = ChainMeans
    
    # MaxVar=numeric(dims)
    # for (i in 1:(MaxMig+1))
    # {
    #   for(j in 1:dims)
    #   {
    #     if(is.na(var(ChainMeans[j,i,],na.rm=TRUE))==FALSE)
    #     {
    #       if(var(ChainMeans[j,i,],na.rm=TRUE)>MaxVar[j])
    #       {
    #         MaxVar[j]=var(ChainMeans[j,i,],na.rm=TRUE)
    #       }
    #     }
    #   }
    # }
    
    #  for(j in 1:dims)
    #      {
    #          if(sqrt(MaxVar[j])>0.05*(max(Coordinates[,j])-min(Coordinates[,j])))
    #             {
    #                 writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations")
    #                 break
    #             }
    #  
    #     }
    
    # if(sqrt(MaxVar)>0.1*max(max(Coords[,1])-min(Coords[,1]),max(Coords[,2])-min(Coords[,2])))
    # {
    #   writeLines("NO CLUSTER CONVERGENCE: You need to re-run the sampler with more iterations");
    # }
    # print(MCMCout$RootProbsR)
    
    rootprobs1 = MCMCout$RootProbsR[1:MCMCout$countR]
    rootprobs2 = MCMCout$RootProbsR[(MCMCout$countR+1):(MCMCout$countR*2)]
    # print(rootprobs1/sum(rootprobs1))
    # print(rootprobs2/sum(rootprobs2))
    
    if(max(abs(rootprobs1/sum(rootprobs1)-rootprobs2/sum(rootprobs2)))>0.5/MCMCout$countR)
    {
      writeLines("NO ROOT CONVERGENCE: You need to re-run the sampler with more iterations");
    }
    
    MCMCout$RootProbsR=MCMCout$RootProbsR[1:MCMCout$countR]+MCMCout$RootProbsR[(MCMCout$countR+1):(MCMCout$countR*2)]
    MCMCout$TreeEdgesR=BPEC.TreeEdges(MCMCout)
    
    SeqLabels = MCMCout$SeqsFileR[MCMCout$SeqLabelsR]
    if(which.max(MCMCout$RootProbsR) <= length(SeqLabels))
    {
      root = SeqLabels[which.max(MCMCout$RootProbsR)]
      writeLines(paste("The most likely root node is ",root,sep=""))
    }
    if(which.max(MCMCout$RootProbsR)>length(SeqLabels))
    {
      root = which.max(MCMCout$RootProbsR)
      writeLines(paste("The most likely root node is extinct",sep=""))
    }
    MCMCout$RootProbsR = MCMCout$RootProbsR/sum(MCMCout$RootProbsR)
    
    UniqueLocs = Coordinates[!duplicated(Coordinates[1:2]),]
    RootProbsUnique = MCMCout$RootLocProbsR
    MCMCout$LocNoR = length(UniqueLocs)
    
    for (i in 1:length(UniqueLocs[,1]))
    {          
      IdenticalLocs = (Coordinates[,1]==UniqueLocs[i,1]) & (Coordinates[,2]==UniqueLocs[i,2])
      TempRoot = sum(RootProbsUnique[IdenticalLocs])
      RootProbsUnique[IdenticalLocs] = TempRoot
      if(sum(IdenticalLocs)>1)
      {
        firstLoc = min(which(IdenticalLocs==TRUE))
        IdenticalLocs[firstLoc] = FALSE
        RootProbsUnique[IdenticalLocs] = 0
      }
      
    }
    
    MCMCout$RootLocProbsR = RootProbsUnique
    
    Rootlocs = numeric(3)
    Rootlocs[1] = sapply(sort(MCMCout$RootLocProbs, index.return=TRUE), `[`, length(MCMCout$RootLocProbs)-1+1)[2]
    Rootlocs[2] = sapply(sort(MCMCout$RootLocProbs, index.return=TRUE), `[`, length(MCMCout$RootLocProbs)-2+1)[2]
    Rootlocs[3] = sapply(sort(MCMCout$RootLocProbs, index.return=TRUE), `[`, length(MCMCout$RootLocProbs)-3+1)[2]
    writeLines(paste("The most likely ancestral locations are ",Rootlocs[1],",",Rootlocs[2],",",Rootlocs[3],sep=""))
    
    HaploIndex = numeric(sum(MCMCout$NoSamplesR))
    counter = 1
    for(i in 1:length(MCMCout$NoSamplesR))
    {
      if(MCMCout$NoSamplesR[i]>0)
      {
        HaploIndex[counter:(counter+MCMCout$NoSamplesR[i]-1)]=i
        counter=counter+MCMCout$NoSamplesR[i]
      }
    }
    
    codatemp = array(0,dim=c(dim(MCMCout$ClusterProbsR)[1],dim(MCMCout$SampleIndicesR)[2]))
    codatemp = MCMCout$ClusterProbsR
    for(i in 1:dim(MCMCout$ClusterProbsR)[1])
    {
      for(j in 1:dim(MCMCout$ClusterProbsR)[2])
      {
        codatemp[i,j]=sum(MCMCout$SampleIndicesR[HaploIndex==i,]==j)
      }
      
      codatemp[i,]=codatemp[i,]/sum(codatemp[i,])
    }
    
    MCMCout$SeqR = MCMCout$SeqR[1:(MCMCout$SeqLengthR*length(RawSeqs))]
    MCMCout$SeqR = array(MCMCout$SeqR,dim=c(length(RawSeqs),MCMCout$SeqLengthR))
    
    MCMCout$ClusterProbsR=codatemp
    # print(MCMCout$ClusterProbsR)
    
    #  if(coda == 1)
    #  {
    dim1 = length(MCMCout$SampleClusterCodaR[,1,1])
    dim2 = length(MCMCout$SampleClusterCodaR[1,,1])
    dim3 = length(MCMCout$SampleClusterCodaR[1,1,])
    chain1 = array(MCMCout$SampleClusterCodaR[,,1:(dim3/2)],dim=c(dim1*dim2,dim3/2))
   
    chain1 = rbind(chain1,MCMCout$SampleRootCodaR[1:(dim3/2)])
    chain2 = array(MCMCout$SampleClusterCodaR[,,(dim3/2+1):(dim3)],dim=c(dim1*dim2,dim3/2))
    chain2 = rbind(chain2,MCMCout$SampleRootCodaR[(dim3/2+1):(dim3)])
    MCMCout$CodaInput=list()
    MCMCout$CodaInput$line1 = mcmc(t(chain1), start = 1, end = dim3/2, thin = 1)
    MCMCout$CodaInput$line2 = mcmc(t(chain2), start = 1, end = dim3/2, thin = 1)
    MCMCout$CodaInput = mcmc.list(MCMCout$CodaInput)
    # }
    
    MCMCout$SampleClusterCodaR = NULL
    MCMCout$ModeInitial = NULL
    MCMCout$ancestralR = NULL
    MCMCout$maxmigR = NULL
    MCMCout$codaR = NULL
    MCMCout$PostSamplesR = NULL
    MCMCout$iterR = NULL
    MCMCout$dsR = NULL
    MCMCout$NSeqR = NULL
    MCMCout$MaxLocR = NULL
    MCMCout$seedsR = NULL
    MCMCout$rootR = NULL
    MCMCout$rootlocsR = NULL
    MCMCout$clustersR = NULL
    MCMCout$SeqCorrR = NULL
  }
  return(MCMCout)
}

