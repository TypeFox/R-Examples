Simpsons <-
function(X,Y,clusterid,clustervars,data,nreps=5000)
{
  if (missing(X)|missing(Y)) stop("'X' and 'Y' must be specified")
  #  if (!require("mclust")) stop("'mclust' package missing")
  #  installs required package  if(!("mclust" %in% .packages(all.available=TRUE))) install.packages("mclust")
  
  
  # Logical to test whether clusterid has been supplied by user:
  EstCluster <- missing(clusterid)
  
  # Logical to test whether clustervars is assigned:
  ClustGiven <- !missing(clustervars)
  
  # collects argumentnames:
  mc <- match.call()
  namex=as.character(mc[['X']]) #defining the X label for the plot
  namey=as.character(mc[['Y']]) # defining the Y label for the plot
  mclustanalysis=NULL
  # If data has been specified
  if (!missing(data))
  {
    # Check whether it is a dataframe
    if (!is.data.frame(data)) stop("'data' not a data frame")
    ## If the data is not a data frame, use the command as.data.frame(data)        
    # Check whether variables exist in dataframe
    ## ik heb hier namey=as.character(mc[['Y']]) gebruikt, anders voert ie de functie steeds uit
    if (is.null(data[[namex]])) stop(paste("Variable",namex,"does not exist in dataframe"))
    if (is.null(data[[namey]])) stop(paste("Variable",namey,"does not exist in dataframe"))
    
    alldata <- data.frame(
      X = data[[namex]],
      Y = data[[namey]])
    
    # Clusterid:
    if (!EstCluster)
    {
      if (is.null(data[[as.character(mc[['clusterid']])]])) stop(paste("Variable",as.character(mc[['clusterid']]),"does not exist in dataframe"))
      
      alldata$clusterid <- data[[as.character(mc[['clusterid']])]]
    }
  } else 
  {
    alldata <- data.frame(X = X, Y = Y)
    
    if (!EstCluster)
    {
      alldata$clusterid <- clusterid
    }
  }
  
  # Estimate clusterID:
  
  if (EstCluster)
  {
    # Cluster variables given:
    if (ClustGiven)
    {
      # If data are specified:
      if (!missing(data))
      {
        varNames <- clustervars
        
        # Check whether variable names have been specified:
        for (var in varNames)
        {
          if (is.null(data[[var]])) stop(paste("Variable",var,"does not exist in dataframe"))
        }
        
        ClustData <- data[varNames]
        
      } else #
      {
        ClustData <- as.data.frame(do.call("cbind",clustervars))
      }
      mclustanalysis <- suppressWarnings(Mclust(ClustData))
    } else
    {
      mclustanalysis <- suppressWarnings(Mclust(alldata[c("X","Y")]))
    }
    
    Nclusters <- mclustanalysis$G            
    alldata$clusterid <- as.vector(mclustanalysis$classification)    
  }
  
  Allclusters<-list()
  Nclusters=max(alldata$clusterid)
  Npercluster<-matrix()
  Allint<-matrix()
  Allbeta<-matrix()
  pvalues<-matrix()
  
  #The following section runs a regression within each cluster, and saves the intercepts and betas
  for(i in 1:Nclusters)
  {
    x=subset(alldata,clusterid==i)
    assign(paste("C", i, sep = ""),x)
    Allclusters[[i]]=x
    Npercluster[i]=nrow(x)
    reg=lm(x[,2]~x[,1])
    Allint[i]=reg$coefficients[1]
    Allbeta[i]=reg$coefficients[2]
	pvals=summary(reg)
	pvals=pvals$coefficients
	pvalues[i]=pvals[8]
  }
  
  #plotting all clusters and drawing regression lines
  plot(alldata[,1],alldata[,2],col = alldata[,3]+1,pch=alldata[,3],xlab=namex,ylab=namey)
  Allbeta=as.numeric(Allbeta)
  Allint=as.numeric(Allint)
  groupreg=lm(alldata[,2]~alldata[,1])
  groupint=groupreg$coefficients[1]
  groupbeta=groupreg$coefficients[2]
  pvalsgr=summary(groupreg)
  pvalsgr=pvalsgr$coefficients
  pvaluesgr=pvalsgr[8]
  for(j in 1:Nclusters)
  {
    a=Allint[j]
    b=Allbeta[j]
    abline(a,b,lty=j,col=j+1)
  }
  abline(groupint,groupbeta,lwd=3,col=1)
  xpd <- par("xpd")
  par(xpd=TRUE)
  legend(mean(par("usr")[c(1,2)]),par("usr")[4],paste("cluster", 1:Nclusters),col= 1:Nclusters+1,pch=1:Nclusters,bty="n",horiz=TRUE,xjust=0.5,yjust=0,cex=1)
  par(xpd=xpd)
  
  #permutation testing of null hypothesis
  totaln=nrow(alldata)  #total n
  permutationtest=matrix(NA,nrow=nreps,ncol=Nclusters) #matrix which will contain every beta estimate for every permuted cluster regression
  lowerboundary=matrix(NA,ncol=Nclusters) #matrix which will contain the lower cut off for each cluster for a significant deviation from the permuted null hypothesis
  upperboundary=matrix(NA,ncol=Nclusters)  #matrix which will contain the lower cut off for each cluster for a significant deviation from the permuted null hypothesis
  
  if (Nclusters>1) message(paste(Nclusters,"clusters detected"))
  if (Nclusters==1) message(paste(Nclusters,"cluster detected"))
  
  for (i in 1:Nclusters)
  {
    clustdata=subset(alldata, clusterid==i) # select ith cluster
    clustern=nrow(clustdata) #define size of ith cluster
    if (Nclusters > 1)
    {
      message(paste("Permuting cluster",i))
      pb <- txtProgressBar(max=nreps,style=3)
      for (j in 1:nreps)
      {
        permid=rep(0:1,c(totaln-clustern,clustern)) 
        permid=sample(permid,size=totaln,replace=F)	#permute cluster id's: Randomly assign a 1 to an equal number of cases as the size of cluster i
        permalldata=cbind(alldata,permid) #append new clusterid to data matrix
        permclust=subset(permalldata,permid==1) #choose subset of data as permuted cluster
        lmclusterbeta=lm(permclust[,2]~permclust[,1]) #run regression within this permuted cluster
        clustbeta=lmclusterbeta$coefficients[2]	 #extract regression beta for 1 permutation of this cluster i
        lmgroupbeta=lm(permalldata[,2]~permalldata[,1]) #run regression for whole dataset (same every iteration but 			defined within loop) 
        groupbeta=lmgroupbeta$coefficients[2] #extract regression beta for whole dataset
        nulldif = clustbeta-groupbeta	#extract difference between group regression and cluster regression
        permutationtest[j,i]= nulldif #write away permuted betas
        setTxtProgressBar(pb, j)
      }						
      close(pb)
        
      lowerboundary[i]=quantile(permutationtest[,i],0.025)
      upperboundary[i]=quantile(permutationtest[,i],0.975)
    }
  }     
  
  
  #looking & reporting for simpsons paradox
  groupsign=sign(groupbeta)
  clustersign=sign(Allbeta)
  
  for (i in 1:Nclusters)
  {  
    group=groupbeta
    clustbeta=Allbeta[i]
    significance=((group-clustbeta)<lowerboundary[i]||(group-clustbeta)>upperboundary[i])
    csign=clustersign[i] #what is the sign of the relationship within the cluster 
    
    if (!is.na(significance))
    {
      if (groupsign!=csign&significance) 
      {
        message(paste("Sign reversal: Simpson's Paradox! Cluster ",i," is significantly different and in the opposite direction compared to the group!",sep=" "))
      } else if (significance) 
      {
        message (paste("Warning: Beta regression estimate in Cluster ",i," is significantly different compared to the group!",sep=""))  
      } else {
        message (paste("No evidence for Simpson's Paradox: Cluster",i,"is not significantly different compared to the group regression",sep=" "))
      }
    } else {
      message("One cluster detected: No evidence for Simpson's Paradox")
    }
    
  }
  


  #   ls()
  results=list(
    Nclusters = Nclusters, 
    clustersize = tapply(alldata$clusterid,alldata$clusterid,length),
    alldata = alldata,
    Allbeta = Allbeta,
    Allint = Allint,
    permutationtest = permutationtest,
    namex = namex,
    namey = namey,
    pvalues=pvalues,
    groupbeta = groupbeta,
    groupint = groupint,
    groupint = groupint,
    pvaluesgr = pvaluesgr,
    totaln = totaln,
    mclustanalysis=mclustanalysis
  ) 
## Here the lists elements get named 
  class(results)<- "Simpson"
  return(results)



}



