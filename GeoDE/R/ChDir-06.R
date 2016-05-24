##########################################
#
# A wrapper function for performing characteristic
# direction analysis
#
############################################

gmt <- NULL
AllGMTfiles <- NULL

chdirAnalysis <- function(datain, sampleclass, gammas=list(1.0), nnull=10,CalculateSig=FALSE)
{  # Test that the inputs are in the correct form and are self-consisent
  
  if(length(sampleclass)!=(length(datain)-1)) stop("number of elements in sampleclass is inconsistent with input data")
  if(!is.data.frame(datain)) stop("Input data is not in the form of a data frame")
  if(FALSE%in%(c("1","2")%in%levels(sampleclass))) stop ("sample class does not include \'1\' and \'2\'")
  if(length(datain[sampleclass==1])<2) stop ("too few controll samples")
  if(length(datain[sampleclass==2])<2) stop ("too few samples")
  
  # Calculate the characteristic direction
  chdirresults <- chdirSig(datain,sampleclass,gammas,nnull=nnull,CalculateSig=CalculateSig)
  
  # produce plots
#  print("plotting...")
  chdirplots(chdirresults,sampleclass,gammas,CalculateSig)
#  print("generating output...")
  
  # Generate result dataframe
  outAll <- lapply(chdirresults[[1]], function(x) {x[sort.list(x^2,decreasing=TRUE),]})
  if(CalculateSig)
  {
    outSig <- mapply( function(x,ns) {
      
      x[sort.list(x^2,decreasing=TRUE)[1:ns],]}, chdirresults[[1]],chdirresults[[6]],SIMPLIFY=FALSE)
    list(chdirprops=chdirresults,results=outSig)  
  }else{list(chdirprops=chdirresults,results=outAll)}
  
}






##########################################
#
# A function to display he results of the characteristic
# Direction function
#
##########################################


chdirplots <- function(chdirresults,sampleclass,gammas=list(1.0),CalculateSig=FALSE)
{
  # Plot the 2-d projections of the data and the characteristic directions
  
  # Plot 2D projections
  mapply(function(a,g)
  {
    arrowscale=0.5
    
    par(las = 1, mar=c(5,4,4,2)+0.1)
    
    plot(chdirresults[[2]][,1],chdirresults[[2]][,2], pch=as.numeric(sampleclass),xlab="PC1", ylab="PC2",main=paste("2D projection of CD: gamma = ", as.character(g)))
    arrowstart <- colMeans(chdirresults[[2]][sampleclass==1,])
    arrows(x0=arrowstart[[1]],y0=arrowstart[[2]],x1=arrowstart[[1]]+arrowscale*a[[1]], y1=arrowstart[[2]]+arrowscale*a[[2]])
  },chdirresults[[3]],gammas)
  
  
  # Plot coeffs of top genes
  
  
  # Plot significance curve
  
  if(CalculateSig)
  {
    mapply(function(a,g)
    {
      plot(a, type="l",xlab="rank",ylab="s", main=paste("Significance curve, gamma = ",g))
    },chdirresults[[5]],gammas)
  }

  
  # Bar plot of chdir components
  
  # Plot 2D projections
  mapply(function(a,g)
  {
        
    par(las = 2, mar=c(7,5,4,2)+0.1)
    

    
    nbars <- min(40,length(a))
    
    
    barplot(t(a[sort.list(a^2,decreasing=TRUE)[1:nbars],]),main=paste("Top genes: gamma = ", as.character(g)),ylab="coeficient")

  },chdirresults[[1]],gammas)
  
}


####################################
# Characteristic direction function WITH significance test #
# Input:
# 1) expression data with genes down the rows and samplees across the columns
# 2) a factors describing which samples are controll (1) and which are perturbed (2)
# 3) a list of gamma values (default (1.0))
# Output
# 1) The chdir in list form
####################################

chdirSig <- function(data,sampleclass,gammas=list(1.0),nnull=10,CalculateSig=FALSE)
{
  
  #
  # Calculate the pca
  #
  
  pca1 <- prcomp(t(as.matrix(data[-1])))
  
  meanvec <- rowMeans(as.matrix(data[-1][sampleclass==2]))-rowMeans(as.matrix(data[-1][sampleclass==1]))
  
  n1 <- sum(sampleclass==1)
  n2 <- sum(sampleclass==2)
  
  #
  #   In order to prevent singularities later we now select the number 
  #of PCs to keep
  #
  
  cumsum <- pca1$sdev^2/sum(pca1$sdev^2)
  keepPC <- length(cumsum[cumsum>0.001])
  
  #
  # Now calculte the characteristic direction in pca space
  #
  V <- pca1$rotation[,1:keepPC]
  R <- pca1$x[,1:keepPC]
  
  Dd <- (t(R[sampleclass==1,])%*%R[sampleclass==1,]+t(R[sampleclass==2,])%*%R[sampleclass==2,])/(n1+n2-2)
  
  sigma <- mean(diag(Dd))
  
  ShrunkMats <- lapply(gammas, function(x) solve(x*Dd + sigma*(1-x)*diag(keepPC)))
  
  # Map back to the full expression space
  
  
  b <- lapply(ShrunkMats, function(x) matrix(V%*%x%*%t(V)%*%meanvec,dimnames=list(c(as.list(as.character(data[[1]]))), 1)))
  
  # Normalize the characteristic directions
  
  
  b <- lapply(b, function(x) x/sqrt(sum(x^2)))
  
  # b<-as.vector(b[[1]])  
  #The projection of the characteristic direction in the first two PCs
  
  
  b2dscale <- colMeans(R[sampleclass==2,1:2])- colMeans(R[sampleclass==1,1:2])
  
  
  b2dscale <- sqrt(sum(b2dscale^2))
  
  #print(b2dscale)
  
  projchdir2d <-lapply(b, function(x) list( b2dscale*as.numeric(as.vector(x)%*%as.vector(V[,1])), b2dscale*as.numeric(as.vector(x)%*%as.vector(V[,2]))))  
  
  
 
  if (CalculateSig) 
  {
    
    
    ########################################
    # Generate a null distribution of chdirs
    ########################################
   
    
    
    
    nu<-n1+n2-2
    y1 <- t(t(mvrnorm(nnull, rep(0, as.numeric(keepPC)), Dd) *sqrt(nu / rchisq(nnull, nu))))
    y2 <- t(t(mvrnorm(nnull, rep(0, as.numeric(keepPC)), Dd) *sqrt(nu / rchisq(nnull, nu))))
    bmeanvec <- colMeans(R[sampleclass==2,])- colMeans(R[sampleclass==1,])
    y <- t(y2-y1)
    
    
    #  y <- mvrnorm(nnull, rep(0, as.numeric(keepPC)), Dd)
    #  y2<-sqrt(nu / rchisq(nnull, nu))
    
    ###############################################
    # For each value of gamma and each of the null
    ###############################################  
#    print("1/2")
    pb <- txtProgressBar(min = 0, max = nnull,style=3) 
    nullchdirs <- lapply(gammas, function(x)
    {
      rowMeans( 
        as.data.frame(
          mapply( function(mv,count)
          {
            
        
            setTxtProgressBar(pb, count)
          
            
            sm <- solve(x*Dd + sigma*(1-x)*diag(keepPC))        
            bn <-  matrix(V%*%sm%*%t(V)%*%as.numeric(mv%*%t(V)),dimnames=list(c(as.list(as.character(data[[1]]))), 1))
            bn <- bn/sqrt(sum(bn^2))
            bn<-bn^2
            bn<-sort(bn,decreasing=TRUE)
          
           
          },as.data.frame(y),c(1:length(y))
          )
        )
      )
    }
    )
    close(pb)  
  

#  print(class(nullchdirs))
#  print(length(nullchdirs))
#  print(dim(nullchdirs))

    ###########################################
    # The ratio to null
    ###########################################
#    print("2/2")
#    pb <- txtProgressBar(min = 0, max = nnull,style=3)
    
    
    ratio <- mapply(function(ba,bn,count) {
#      setTxtProgressBar(pb, count)
      relerr <- sort(ba^2,decreasing=TRUE)/bn
      
      relerr <- cumsum(relerr)/sum(relerr)-c(1:length(meanvec))/length(meanvec)
    }, b,nullchdirs, c(1:length(nullchdirs)),SIMPLIFY=FALSE)
#    close(pb)
#    print("This may take several minutes, please wait...")
    nsiggenes <-lapply(ratio,function(x) which.max(x))
      
    list(chdir=b,pca2d=R[,1:2],chdir_pca2d=projchdir2d,null_rank_dist=nullchdirs,ratio_to_null=ratio,number_sig_genes=nsiggenes)
    
    
  }else
  {
    list(chdir=b,pca2d=R[,1:2],chdir_pca2d=projchdir2d) 
  } 
  
  
}






#################################
#
#
# Multi GMT PAEA analysis
#
###################################

multigmtPAEAAnalysis <- function(chdirresults,gmtfiles=AllGMTfiles,gammas=c(1.0),casesensitive=FALSE,showprogress=TRUE)
{

  pb <- txtProgressBar(min = 0, max = length(gmtfiles),style=3)
  
  mapply(function(gmtfile,count)
    {
    print(count)
    setTxtProgressBar(pb, count)

    data(list=gmtfile, envir = environment())
    

    ar <- PAEAAnalysis(chdirresults,gmt,gammas=gammas,casesensitive=casesensitive,showprogress=FALSE)
    
    write.table(ar$p_values, paste("PAEA_p_values-",gmtfile,".txt",sep=""), sep="\t")
    
  }, gmtfiles,c(0:(length(gmtfiles)-1)))

  close(pb)
  
}



#################################################
#
# A wrapper function for PAEA
#
#################################################

PAEAAnalysis <- function(chdirresults,gmtfile,gammas=c(1.0),casesensitive=FALSE,showprogress=TRUE)
{
  # Extract the names of the gmt lines and the gmt lines
  gmtlinenames <- lapply(gmtfile, function(x) x[[1]])
  gmtlines <- lapply(gmtfile, function(x) x[-1])
  
#  print("processed gmt")
  
  # Calculate the PEAE results for each gmt line
  if(showprogress){
    pb <- txtProgressBar(min = 0, max = length(gmtlines),style=3)
    PAEAresults <-mapply(function(x,count) 
      {
      setTxtProgressBar(pb, count)
      PAEA(chdirresults[[1]],x,casesensitive=casesensitive)
    },gmtlines,c(1:length(gmtlines)),SIMPLIFY=FALSE)    
    close(pb)
    
  }else{
   
    PAEAresults <-lapply(gmtlines, function(x) PAEA(chdirresults[[1]],x,casesensitive=casesensitive))
    
  }
  


  
  gammalabels <-unlist(lapply(gammas, function(x) paste("gamma=",x)))
  
  # Sort the lines and output the results
  
  pvalues<-lapply(PAEAresults, function(x) x[[2]])

  

  pvalues<-matrix(unlist(pvalues),ncol=length(gmtlines), dimnames=list(gammalabels,gmtlinenames))
  
  pavalues<-lapply(PAEAresults, function(x) x[[1]])
  pavalues<-matrix(unlist(pvalues),ncol=length(gmtlines), dimnames=list(gammalabels,gmtlinenames))
  

#  print(class(pvalues))
#  print(length(pvalues))
#  print(dim(pvalues))
  
  # identify the ordering by the first gamma value
  gmtp <- sort.list(pvalues[1,])
  
  # Produce a bar-graph plot of the top enriched lists
  nbars<-min(10,length(gmtlines))
  
  
#  print(pvalues[1,gmtp[1:nbars]])
 
  pvalues <- pvalues+1e-200

#  print(pvalues[1,gmtp[1:nbars]])

  par(las = 1, mar=c(5,20,1,1))
  barplot(-rev(log(pvalues[1,gmtp[1:nbars]])), horiz = TRUE,xlab="-log(p)")
  
  list(p_values=t(pvalues[,gmtp]),principal_angles=t(pavalues[,gmtp]))
}


#####################################
# PAEA function
# Input:
# 1) A characteristic direction
# 2) a factor of gene names corrsponding to the elements of the characteristic direction
# 3) a line from a gmt file (class character)
# Output
# 1) The principal angle
# 2) The p value
#####################################
PAEA<-function(chdir,gmtline,casesensitive=FALSE)
{
  genenames<-rownames(chdir[[1]])
  
#  print("in PAEA")

  
  keepgenes <-is.element(toupper(genenames), toupper(gmtline))
  ngenes <-length(genenames)
  
  # A matrix with a 1 in rows corresponding to genes in the gmt line
  #gvec <- as.matrix(lapply(keepgenes, function(x) if (x==TRUE) 1 else 0))
  
  gpos <-which(toupper(genenames)%in%toupper(gmtline))
 
  nset <-length(gmtline)

  n <-ngenes
  m <-length(gpos)

#  print("n genes and set = ")
#  print(n)
#  print(m)


#  print(gpos)
#  print(length(gpos))
#  print("doing gsa")
  
  if((length(gpos)>0)&(m>1)&(m<n))
  {
  gsa <- as.matrix(sparseMatrix(gpos, c(1:length(gpos)), dims=c(ngenes,length(gpos)),x=1))
  
  #   Calculate the principal angle
  
  #  prod <- lapply(chdir, function(x) t(x)%*%gsa)
  
#  print("calculating pa")
  
  principalangle <- lapply(chdir, function(x) svd((t(x)/sqrt(sum(x^2)))%*%gsa, nu=0, nv=0))
  
 
  
  #####
  # Isotropic null distribution
  #####
  
  #2 (Gamma[(n)/2]/(Gamma[(n - m)/2]*Gamma[m/2])) (Sin[theta])^(n - m - 1)*(Cos[theta])^(m - 1);
  

  pac <- function(theta) {2*(1./sqrt(2*pi))*exp((n/2)*log(n/(n - m))+(m/2)*log((n - m)/m)+(1/2)*log(m/(2*n)*(n - m))+(n-m-1)*log(sin(theta))+(m-1)*log(cos(theta)))}
  
  
  #  Ncdfpa <-integrate(pac,0,as.numeric(principalangle[[1]])) 
  #  list(as.numeric(principalangle[[1]]),Ncdfpa)
  
  Ncdfpa <-lapply(principalangle,function(pa) integrate(pac,0,acos(as.numeric(pa)))$value)
  
  list(principal_angle=principalangle,p_val=Ncdfpa)
  
  }else
  {
    principalangle <- lapply(chdir, function(x) 0.) 
    Ncdfpa <-lapply(principalangle,function(pa) 1.0)
    
    list(principal_angle=principalangle,p_val=Ncdfpa)
  }
  
  
}
