
# Generating data sets via a factorial design (Qiu and Joe, 2006), 
# which has factors degree of separation, number of clusters, 
# number of non-noisy variables, number of noisy variables, 
# number of outliers, number of replicates, etc. 
# The separation between any cluster and its nearest neighboring clusters 
# can be set to a specified value. The covariance matrices of clusters 
# can have arbitrary diameters, shapes and orientations.

#   Qiu, W.-L. and Joe, H. (2006)
#   Generation of Random Clusters with Specified Degree of Separaion.
#   \emph{Journal of Classification}, \bold{23}(2),  315--334.
#
#  Factors                          Levels
#---------------------------------------------
# Number of clusters                3, 6, 9
# Degree of separation              close, separated, well-separated
# Number of non-noisy variables     4, 8, 20
# Number of noisy variables         1, 0.5p1, p1
#---------------------------------------------
# p1 is the number of non-noisy variables
#
# sepVal -- Levels of the degree of separation. 
#          It has been "calibrated" against two univariate 
#          clusters generated from N(0, 1) and N(0, A), respectively,  
#          With A=4, sepVal=0.01 indicates a close cluster structure. 
#          With A=6, sepVal=0.21 indicates a separate cluster structure. 
#          With A=8, sepVal=0.34 indicates a well-separated cluster structure.
#      'close' -- 'sepVal=0.01'
#      'separated' -- 'sepVal=0.21'
#      'well-separated' -- 'sepVal=0.342'
# sepLabels -- labels for 'close', 'separated', and 'well-separated' 
# numClust -- levels of numbers of clusters
# numNonNoisy -- levels of numbers of non-noisy variables
# numNoisy -- levels of numbers of noisy variables
# numOutlier -- number/ratio of outliers. If numOutlier is a positive integer,
#   then numOutlier means the number of outliers. If numOutlier is a real
#   number in (0,1), then numOutlier means the ratio of outliers, i.e. 
#   the number of outliers is equal to round(numOutlier*N), 
#   where N is the total number of non-outliers.
# numReplicate -- the number of data sets to generate for each combination.
#   the default value is 3 as in the design in Milligan (1985) 
#   (An algorithm for generating artificial test clusters,
#   Psychometrika, 50, 123-127)
# fileName -- the output data file names have the format 
#       [fileName]J[j]G[g]v[p1]nv[p2]out[numOutlier]_[numReplicate].[extension] 
#   where
#      'extension' can be 'dat', 'log', 'mem', or 'noisy',
#      'J' means separation index, 'G' means number of clusters,
#      'v' means the number of non-noisy variables,
#      'nv' means the number of noisy variables,
#      'out' means the number of outliers
# clustszind -- cluster size indicator. 
#      clustszind=1 indicates that all cluster have equal size.  
#         The size is specified by clustSizeEq.
#      clustszind=2 indicates that the cluster sizes are randomly 
#         generated from the range rangeN.
#      clustszind=3 indicates that the cluster sizes are specified
#         via the vector clustSizes.
#      The default value is 2 so that the generated clusters are more realistic.
# clustSizeEq -- if clustszind=1, then this is the constant cluster size 
#      The default value 100 is a reasonable cluster size.
# rangeN -- if clustszind=2, then the cluster sizes are randomly generaged
#         from the range 'rangeN'. The default range is [50, 200] which
#         can produce reasonable variability of cluster sizes.
# clustSizes -- if clustszind=3, then the cluster sizes are specified via
#         clustSizes. An input is required with
#         the default value of 'clustSizes' set as 'NULL'.
# covMethod -- methods to generate random covariance matrices. 
#      'eigen' method first generates the eigenvalues of the covariance matrix,
#          then generate eigenvectors to construct the covariance matrix. 
#      'unifcorrmat' method first generates a correlation 
#          matrix via the method proposed by 
#           Joe H (2006). Generating random correlation matrices based on 
#           partial correlations. J. Mult. Anal. Vol. 97, 2177--2189
#          Then, it generate variances from the range 'rangeVar' to
#          construct the covariance matrix. 
#
#      'onion' method
# extension of onion method, using Cholesky instead of msqrt
#Ghosh, S., Henderson, S. G. (2003). Behavior of the NORTA method for
#correlated random vector generation as the dimension increases.
#ACM Transactions on Modeling and Computer Simulation (TOMACS)
#v 13 (3), 276-294.
#
#      'c-vine' method
#      # Random correlation with the C-vine (different order of partial
# correlations). Reference for vines: Kurowicka and Cooke, 2006,
# Uncertainty Analysis with High Dimensional Dependence Modelling,
# Wiley, 2006.
#
#      The default method is 'eigen' so that the user can directly 
#      specify the range of the 'diameters' of clusters.
# rangeVar -- if 'covMethod="unifcorrmat"', then to generate a covariance 
#      matrix, we first generate a correlation matrix, then randomly generate
#      variances from the range specified by 'rangeVar'.
#      The default range is [1, 10] which can generate reasonable
#      variability of variances.
# lambdaLow -- if 'covMethod="eigen"', when generating a covariance matrix, 
#      we first generate its eigenvalues. The eigenvalues are randomly 
#      generated from the interval [lambdaLow, lambdaLow*ratioLambda]. 
#      Note that lambdaLow should be positive.
#      In our experience, the range [lambdaLow=1, ratioLambda=10]
#      can give reasonable variability of the diameters of clusters.
# ratioLambda -- if 'covMethod="eigen"', lambdaUpp=lambdaLow*ratioLambda, 
#      is the upper bound of the eigenvalues of a covariance matrix
#      and lambdaLow is the lower bound.
#      In our experience, the range [lambdaLow=1, lambdaUpp=10]
#      can give reasonable variability of the diameters of clusters.
# alphad -- parameter for c-vine and onion method to generate random correlation matrix
#       eta=1 for uniform. eta should be > 0
# eta -- parameter for c-vine and onion method to generate random correlation matrix
#       eta=1 for uniform. eta should be > 0
# rotateind-- if rotateind =TRUE, then rotate data so that we may not detect the
#      full cluster structure from scatterplots, otherwise do not rotate.
#      By default, 'rotateind=TRUE' to generate more realistic data sets.
# iniProjDirMethod -- indicating the method to get initial projection direction 
#      By default, the sample version of the Su and Liu (SL) projection 
#      direction ((n1-1)*S1+(n2-1)*S2)^{-1}(mu2-mu1) is used,
#      where mui and Si are the mean vector and covariance
#      matrix of cluster i. Alternatively, the naive projection
#      direction (mu2-mu1) is used. 
#      Su and Liu (1993) JASA 1993, vol. 88 1350-1355.
# projDirMethod -- takes values "fixedpoint" or "newton"
#      "fixedpoint" means that we get optimal projection direction
#        by solving the equation d J(w) / d w = 0, where
#        J(w) is the separation index with projection direction 'w'
#      "newton" method means that we get optimal projection driection
#        by the method proposed in the appendix of Qiu and Joe (2006) 
#        "Generation of random clusters with specified degree of separation",
#        Journal of Classification Vol 23(2), 315-334.
# alpha -- tuning parameter for separation index to indicating the percentage 
#      of data points to downweight. We set 'alpha=0.05' like we set
#      the significance level in hypothesis testing as 0.05.
# ITMAX -- when calculating the projection direction, we need iterations. 
#      ITMAX gives the maximum iteration allowed.
#      The actual #iterations is usually much less than the default value 20.
# eps -- A small positive number to check if a quantitiy \eqn{q} is 
#      equal to zero. If \eqn{|q|<}\code{eps}, then we regard \eqn{q} 
#      as equal to zero.  \code{eps} is used to check if an algorithm 
#      converges. The default value is \eqn{1.0e-10}.
# quiet -- a flag to switch on/off the outputs of intermediate results.
#      The default value is 'TRUE'.
# outputEmpirical -- indicates if empirical projection directions
#      and separation indices should be output to files.
#      These information usually are useful to check the cluster
#      structures. Hence, by default, 'outputEmpirical=TRUE'.
# outputInfo -- indicates if separation information dataframe should be output.
#      The file name extension is .log
simClustDesign<-function(numClust=c(3,6,9), 
                         sepVal=c(0.01, 0.21, 0.342), 
                         sepLabels=c("L", "M", "H"), 
                         numNonNoisy=c(4,8,20), 
                         numNoisy=NULL, 
                         numOutlier=0, 
                         numReplicate=3, 
                         fileName="test", 
                         clustszind=2, 
                         clustSizeEq=50, 
                         rangeN=c(50,200), 
                         clustSizes=NULL,
                         covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"), 
                         rangeVar=c(1, 10), 
                         lambdaLow=1, 
                         ratioLambda=10, 
                         alphad=1, 
                         eta=1, 
                         rotateind=TRUE, 
                         iniProjDirMethod=c("SL", "naive"), 
                         projDirMethod=c("newton", "fixedpoint"), 
                         alpha=0.05, 
                         ITMAX=20, 
                         eps=1.0e-10, 
                         quiet=TRUE, 
                         outputDatFlag=TRUE,
                         outputLogFlag=TRUE,
                         outputEmpirical=TRUE, 
                         outputInfo=TRUE)
{
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  covMethod<-match.arg(arg=covMethod, choices=c("eigen", "onion", "c-vine", "unifcorrmat"))
  numClust<-as.integer(numClust)

  # checks for valid inputs
  if(prod(numClust<1) || !is.integer(numClust)) 
  { 
    stop("The number 'numClust' of clusters should be a positive integer!\n")
  }
  numNonNoisy<-as.integer(numNonNoisy)
  if(prod(numNonNoisy<2) || !is.integer(numNonNoisy)) 
  { 
    stop("The number 'numNonNoisy' of non-noisy variables should be an integer greater than 1!\n") 
  }
  if(prod(sepVal<= -0.999) || prod(sepVal >= 0.999)) 
  { 
    stop("The desired separation index 'sepVal' should be in the range (-0.999, 0.999)!\n")
  }
  numReplicate<-as.integer(numReplicate)
  if(numReplicate<1 || !is.integer(numReplicate))
  { 
    stop("The number 'numReplicate' should be a positive integer!\n")
  }
  if(!is.null(numNoisy))
  {
    numNoisy<-as.integer(numNoisy)
    if(numNoisy<0 || !is.integer(numNoisy)) 
    { 
      stop("The number 'numNoisy' of noisy variables should be a non-negative integer!\n") 
    }
  }
  if(numOutlier<0)
  {
    stop("'numOutlier' should be positive!\n")
  }
  if(!is.element(clustszind, c(1,2,3))) 
  { 
    stop("Cluster size indicator 'clustszind' should be 1, 2, or 3!\n")
  }
  clustSizeEq<-as.integer(clustSizeEq)
  if(clustSizeEq<2 || !is.integer(clustSizeEq)) 
  { 
    stop("Cluster size 'clustSizeEq' should be an integer greater than 1!\n") 
  }
  rangeN<-as.integer(rangeN)
  if(length(rangeN)!=2) 
  { 
    stop("The range 'rangeN' for cluster sizes should be a numeric vector of length 2!\n") 
  }
  if(rangeN[1]>rangeN[2]) 
  { 
    stop("First element of 'rangeN' should be smaller than second!\n") 
  }
  if(rangeN[1]<2 || !is.integer(rangeN[1])) 
  { 
    stop("The lower bound 'rangeN[1]' for the range of cluster sizes should be an integer greater than 1!\n") 
  }
  if(!is.integer(rangeN[2])) 
  { 
    stop("The upper bound 'rangeN[2]' for the range of cluster sizes should be integer!\n") 
  }
  if(clustszind==3)
  {
    len<-length(clustSizes)
    if(len!=numClust || is.null(clustSizes)) 
    {
      stop("The number of elements in 'clustSizes' should be equal the number of elements in 'numClust' when the value of 'clustszind' is equal to 3!\n")
    }
    clustSizes<-as.integer(clustSizes)
    for(i in 1:len)
    { if(clustSizes[i]<1 || !is.integer(clustSizes[i]))
      { stop(paste("The cluster size for the ", i, "-th cluster should be an integer greater than 1!\n", sep=""))
      }
    }
  }
  if(rangeVar[1]>rangeVar[2]) 
  { 
    stop("First element of 'rangeVar' should be smaller than second!\n") 
  }
  if(rangeVar[1]<0) 
  { 
    stop("The lower bound 'rangeVar[1]' for the range of variances should be positive!\n") 
  }
  if(lambdaLow<0) 
  { 
    stop("The lower bound 'lambdaLow' of eigenvalues of cluster covariance matrices should be greater than zero!\n") 
  }
  if(ratioLambda<1) 
  { 
    stop("The ratio 'ratioLambda' of the upper bound of the eigenvalues to the lower bound of the eigenvalues of cluster covariance matrices should be greater than 1!\n")
  }
  alphad<-as.numeric(alphad)
  if(alphad<0)
  {
    stop("'alphad' should be positive!\n")
  }
  eta<-as.numeric(eta)
  if(eta<0)
  {
    stop("'eta' should be positive!\n")
  }

  if(!is.logical(rotateind))
  {
    stop("The value of the rotation indicator 'rotateind' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<=0 || !is.integer(ITMAX))
  {
    stop("The maximum iteration number allowed 'ITMAX' should be a positive integer!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be between (0, 0.01]!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  if(!is.logical(outputDatFlag))
  {
    stop("The value of the indicator 'outputDatFlag' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  if(!is.logical(outputLogFlag))
  {
    stop("The value of the indicator 'outputLogFlag' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  if(!is.logical(outputEmpirical))
  {
    stop("The value of the indicator 'outputEmpirical' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  if(!is.logical(outputInfo))
  {
    stop("The value of the indicator 'outputInfo' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  # end of checks of valid inputs, loop begins

  loop<-0
  cat("Generating data sets. Please wait ...\n")
   
  for(j in 1:length(sepVal))
  { for(g in numClust)
    { for(p1 in numNonNoisy)
      { if(is.null(numNoisy))
        { numNoisy<-c(1, round(p1/2), p1) } 
        for(p2 in numNoisy)
        { p<-p1+p2 # total number of variables
          loop<-loop+1
          tmpfileName<-paste(fileName,"J", sepLabels[j],"G",g,
            "v", p1, "nv",p2, "out", numOutlier, sep="")
          res<-genRandomClust(numClust=g, numNonNoisy=p1, sepVal=sepVal[j], 
               numNoisy=p2, numReplicate=numReplicate, numOutlier=numOutlier, 
               fileName=tmpfileName,  clustszind=clustszind, 
               clustSizeEq=clustSizeEq, rangeN=rangeN, 
               clustSizes=clustSizes, covMethod=covMethod, 
               rangeVar=rangeVar, lambdaLow=lambdaLow, 
               ratioLambda=ratioLambda,  rotateind=rotateind, 
               iniProjDirMethod=iniProjDirMethod, 
               projDirMethod=projDirMethod, alpha=alpha, 
               ITMAX=ITMAX, eps=eps, quiet=quiet, 
               outputDatFlag=outputDatFlag,
               outputLogFlag=outputLogFlag,
               outputEmpirical=outputEmpirical, 
               outputInfo=FALSE)

          if(loop>1) 
          { infoFrameTheory<-rbind(infoFrameTheory, res$infoFrameTheory)
            if(outputEmpirical)
            { infoFrameData<-rbind(infoFrameData, res$infoFrameData) }
            else { infoFrameData<-NULL }
            datList[[loop]]<-res$datList
            memList[[loop]]<-res$memList
            noisyList[[loop]]<-res$noisyList
          } else {
            infoFrameTheory<-res$infoFrameTheory
            if(outputEmpirical)
            { infoFrameData<-res$infoFrameData }
            else { infoFrameData<-NULL }
           
            datList<-list(res$datList)
            memList<-list(res$memList)
            noisyList<-list(res$noisyList)
          } 
 
        }
      }
    }
  }

  infoFrameTheory<-data.frame(infoFrameTheory)
  nr<-nrow(infoFrameTheory)
  rownames(infoFrameTheory)<-1:nr
  
  if(outputEmpirical)
  { infoFrameData<-data.frame(infoFrameData) 
    nr<-nrow(infoFrameData)
    rownames(infoFrameData)<-1:nr
  }
  else { infoFrameData<-NULL }

  names(datList)<-1:loop
  names(memList)<-1:loop
  names(noisyList)<-1:loop

  if(outputInfo)
  { fileNameInfo<-paste(fileName, "_info.log", sep="")
    msg<-"Theoretical separation information data frame>>>>>>>>>>>\n"
    write.table(msg, append=FALSE, file=fileNameInfo, quote=FALSE,
                row.names=FALSE, col.names=FALSE)
    msg<-colnames(infoFrameTheory)
    write(msg, append=TRUE, file=fileNameInfo, ncolumns=7, sep=" ")
    tt<-infoFrameTheory
    tt[,1:6]<-round(tt[,1:6], 3)
    write.table(tt, file=fileNameInfo, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE, append=TRUE)

    if(outputEmpirical)
    { msg<-"\nEmpirical separation information data frame>>>>>>>>>>>\n"
      write.table(msg, append=TRUE, file=fileNameInfo, quote=FALSE,
                  row.names=FALSE, col.names=FALSE)
      msg<-colnames(infoFrameData)
      write(msg, append=TRUE, file=fileNameInfo, ncolumns=7, sep=" ")
      tt<-infoFrameData
      tt[,1:6]<-round(tt[,1:6], 3)
      write.table(tt, file=fileNameInfo, quote=FALSE, sep="\t",
                  row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }

  cat("The process is completed successfully!\n")
  invisible(list(infoFrameTheory=infoFrameTheory, 
              infoFrameData=infoFrameData, 
              datList=datList, memList=memList, noisyList=noisyList))
} 

