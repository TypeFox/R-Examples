BPEC.CovariatesPlot = function(CovNames,MCMCout,colorcode=c(7,5,6,3,2,8,4,9,10))
{
  writeLines("Creating covariate distribution plot...")
  MeanSamples = MCMCout$SampleMeansR
  CovSamples = MCMCout$SampleCovsR 
  NoClusters=dim(MeanSamples)[2]
  dims=dim(MeanSamples)[1]
  SubSeq=seq(1,length(CovSamples[1,1,1,]),length.out=20)
  MeanSamplesSub = MeanSamples[,,SubSeq]
  
  fullclust=numeric(NoClusters)
  for(i in 1:NoClusters)
  {
    if(length(which(!is.na(MeanSamplesSub[1,i,])))>0)
      #                  if(length(which(!is.na(MCMCout$SampleMeansR[1,i,])))>length(MCMCout$SampleMeansR[1,i,])/2)
    {
      fullclust[i]=1
    }
    else
    {
      fullclust[i]=0
    }
  }
 
  if(dims>2)
  {      
    for(i in 3:dims)
    {
      # w=seq(min(Means[i,])-2*max(sqrt(Covs[i,i,])),max(Means[i,])+2*max(sqrt(Covs[i,i,])),length=100)
      # w=seq(min(MeanSamples[i,,],na.rm=TRUE)-2*max(sqrt(CovSamples[i,i,,]),na.rm=TRUE),max(MeanSamples[i,,],na.rm=TRUE)+2*max(sqrt(CovSamples[i,i,,]),na.rm=TRUE),length=100)
      # w=seq(min(MeanSamples[i,,],na.rm=TRUE)-1*mean(sqrt(CovSamples[i,i,,]),na.rm=TRUE),max(MeanSamples[i,,],na.rm=TRUE)+1*mean(sqrt(CovSamples[i,i,,]),na.rm=TRUE),length=100)
      w=seq(mean(MeanSamples[i,,],na.rm=TRUE)-3*mean(sqrt(CovSamples[i,i,,]),na.rm=TRUE),mean(MeanSamples[i,,],na.rm=TRUE)+3*mean(sqrt(CovSamples[i,i,,]),na.rm=TRUE),length=100)
      maxax=0
      
      for(j in 1:NoClusters)
      {
        if(fullclust[j]==0)
        {
          next
        }                
        for(it in 1:dim(MeanSamples)[3])
        {  
          if(is.na(MeanSamples[i,j,it])==FALSE)
          {
            dnst=dnorm(w,MeanSamples[i,j,it],sqrt(CovSamples[i,i,j,it]))
            if(max(dnst)>maxax)
            {
              maxax=max(dnst)              
            }                 
          }
        }           
      }
      
      plot(w,dnst,col='white',ylim=c(0,maxax),type='l',xlab=CovNames[i],ylab="")     
      
      for(j in 1:NoClusters) 
      {
        if(fullclust[j]==0)
        {
          next
        }
        DensIt=array(NA,dim=c(length(w),dim(MeanSamples)[3]))
        for(it in 1:dim(MeanSamples)[3])
        {                      
          if(is.na(MeanSamples[i,j,it])==FALSE)
          {
            DensIt[,it]=dnorm(w,MeanSamples[i,j,it],sqrt(CovSamples[i,i,j,it]))                        
          }
        }                              
        DensUp=w
        DensLow=w
        DensMed=w
        for(jj in 1:length(w))
        {
          DensUp[jj]=quantile(DensIt[jj,],0.9,na.rm=TRUE)
          DensLow[jj]=quantile(DensIt[jj,],0.1,na.rm=TRUE)
          DensMed[jj]=quantile(DensIt[jj,],0.5,na.rm=TRUE)
        }
        ColCont=col2rgb(colorcode[j], alpha = FALSE)
        polygon(c(w,rev(w)),c(DensLow,rev(DensUp)),col=c(rgb(ColCont[1]/255,ColCont[2]/255,ColCont[3]/255,alpha=0.3)), border = NA)
        lines(w,DensMed,col=colorcode[j])        
      }        
    }
  }
}
