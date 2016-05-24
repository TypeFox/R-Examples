# v1.2.5  
#  (1) fixed a bug in function 'genMemSize'. The 'numClust' inside
#     the function 'genMemSize' should be 'G' 
#
# Generating cluster data sets with specified degree of separation.
# The separation between any cluster and its nearest neighboring cluster 
# can be set to a specified value. 
# The covariance matrices of clusters can have arbitrary diameters, 
# shapes and orientations.

################
# Cluster generating Algorithm.
# (1) The covariance matrices of generated clusters can be arbitrary
#     positive definite matrices.
# (2) The minimum separation index between a cluster and its nearest
#     neighboring cluster is equal to a user-specified value.
#
# References:
#  Joe, H. (2006)
#  Generating Random Correlation Matrices Based on Partial Correlations. 
#  \emph{Journal of Multivariate Analysis}, \bold{97}, 2177--2189.
#
#  Milligan G. W. (1985) 
#  An Algorithm for Generating Artificial Test Clusters.
#  \emph{Psychometrika} \bold{50}, 123--127.
#
#  Qiu, W.-L. and Joe, H. (2006a)
#  Generation of Random Clusters with Specified Degree of Separaion.
#  \emph{Journal of Classification}, \bold{23}(2), 315--334.
#
#  Qiu, W.-L. and Joe, H. (2006b)
#  Separation Index and Partial Membership for Clustering.
#  \emph{Computational Statistics and Data Analysis}, \bold{50}, 585--603.
#
#  Su, J. Q. and Liu, J. S. (1993)
#  Linear Combinations of Multiple Diagnostic Markers.
#  \emph{Journal of the American Statistical Association}, \bold{88}, 1350--1355
#################

#################
# The shifted vertex method is proposed by Dr. Harry Joe. The algorithm is:
#
# Let e1=(1,0,...,0). The first two vertices are -e1 and e1.
# Let the numNonNoisy+1 vertices be labelled v(1),...,v(numNonNoisy+1).
#  numNonNoisy is the number of non-noisy variables.
#  v1=A* (-e1), v2=A*e1
#  vi=sqrt(3)*A*e_{i-1}, i=2,...,numNonNoisy+1
#  sqrt(3) make sure the triangle are equilateral
# One way to get an arbitrary number G of clusters in numNonNoisy dimensions:
# For G<=numNonNoisy+1, use vertices up to v(G).
# For G>numNonNoisy+1, start adding vertices from the following sequence 
#   after v(numNonNoisy+1).
#
# v(2)+2*A*e1, ..., v(numNonNoisy+1)+2*A*e1,
# v(2)+4*A*e1, ..., v(numNonNoisy+1)+4*A*e1,
# v(2)+6*A*e1, ..., v(numNonNoisy+1)+6*A*e1,
# Essentially this just keeps on adding points on a shifted symmetric simplex.
#
# Algorithm: outline of the steps are:
#
#a) input G, numNonNoisy, sepVal, #noisy
#b) generate the G vertices/centers in numNonNoisy dimensions
#c) determine a multiplier A 
#d) generate G covariance matrices
#e) generate a rotation matrix in dimension numNonNoisy
#f) apply the rotation matrix to the G centers and covariance matrices
#g) add noisy variables
#h) apply random permutations to rows and columns of cluster data matrices
#################

#################
# The degree of separation is based on the separation index.
# To calibrate the concepts "well-separated", "separated", and "close",
# we use two clusters from
# two univariate normal distributions N(0, sigma1^2) and N(A, sigma2^2).
# A=4 corresponds to "close" cluster structure;
# A=6 corresponds to "separated" cluster structure;
# A=8 corresponds to "well-separated" cluster structure;
# The corresponding separation indices are:
# 0.01011020, 0.2096862, 0.34229
#################

#################
# The arguments below are also inputs to other functions beside
#   the main function genRandomClust(). The documentation of the 
#   the input arguments will not be repeated.
# numClust -- the number of clusters. 
# numNonNoisy -- the number of non-noisy variables. 
# sepVal -- the minimum separation index specified as a priori. The default 
#          value is 0.01 which is the value of the separation index for
#          two univariate clusters generated from N(0, 1) and N(0, A),
#          respectively, with A=4.
#          With A=4, sepVal=0.01 indicates a close cluster structure. 
#          With A=6, sepVal=0.21 indicates a separate cluster structure. 
#          With A=8, sepVal=0.34 indicates a well-separated cluster structure.
# numReplicate -- the number of data sets to generate for each combination.
#   the default value is 3 as in the design in Milligan (1985) 
#   (An algorithm for generating artificial test clusters,
#   Psychometrika, 50, 123-127)
# numNoisy -- the number of noisy variables. The default values of 'numNoisy'
#   and 'numOutlier' are 0 so that we get 'clean' data sets. 
# numOutlier -- number/ratio of outliers. If numOutlier is a positive integer,
#   then numOutlier means the number of outliers. If numOutlier is a real
#   number in (0, 1), then numOutlier means the ratio of outliers, i.e. 
#   the number of outliers is equal to round(numOutlier*N), 
#   where N is the total number of non-outliers. 
#   If numOutlier is a real number greater than 1, then numOutlier is rounded 
#   to an integer. 
#   The default values of 'numNoisy' and 'numOutlier' are 0 so that we 
#   get 'clean' data sets. 
# fileName --- the fileNames of data sets will start with "fileName" and 
#   followed by numbers, then followed by ".dat".
#   The log, membership, and noisy set files have the same format 
#   except the file extension are ".log", ".mem", and ".noisy" 
#   respectively. The default value is 'test'.
# clustszind -- cluster size indicator. 
#      clustszind=1 indicates that all cluster have equal size.  
#         The size is specified by clustSizeEq.
#      clustszind=2 indicates that the cluster sizes are randomly 
#         generated from the range rangeN.
#      clustszind=3 indicates that the cluster sizes are specified
#         via the vector clustSizes.
#      The default value is 2 so that the generated clusters are more realistic.
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
#        Journal of Classification Vol 23(2), 315--334 
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
#      The file name has the format [fileName]_info.log
#############
genRandomClust<-function(numClust, 
                         sepVal=0.01, 
                         numNonNoisy=2, 
                         numNoisy=0, 
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
  if(numClust<1 || !is.integer(numClust)) 
  { 
    stop("The number 'numClust' of clusters should be a positive integer!\n")
  }
  numNonNoisy<-as.integer(numNonNoisy)
  if(numNonNoisy<2 || !is.integer(numNonNoisy)) 
  { 
    stop("The number 'numNonNoisy' of non-noisy variables should be an integer greater than 1!\n") 
  }
  if(sepVal<= -0.999 || sepVal >= 0.999) 
  { 
    stop("The desired separation index 'sepVal' should be in the range (-0.999, 0.999)!\n")
  }
  numReplicate<-as.integer(numReplicate)
  if(numReplicate<1 || !is.integer(numReplicate))
  { 
    stop("The number 'numReplicate' should be a positive integer!\n")
  }
  numNoisy<-as.integer(numNoisy)
  if(numNoisy<0 || !is.integer(numNoisy)) 
  { 
    stop("The number 'numNoisy' of noisy variables should be a non-negative integer!\n") 
  }
  if(numOutlier<0)
  {
    stop("'numOutlier' should be non-negative!\n")
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
    stop("The range 'rangeN' for cluster size should be a numeric vector of length 2!\n") 
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
    stop("The upper bound 'rangeN[2]' for the range of cluster size should be integer!\n") 
  }
  if(clustszind==3)
  {
    len<-length(clustSizes)
    if(len!=numClust || is.null(clustSizes)) 
    {
      stop("The number of elements in 'clustSizes' should be equal to the number 'numClust' of clusters when the value of 'clustszind' is equal to 3!\n")
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
  # end of checks of valid inputs, algorithm begins


  if(!quiet)
  { 
    cat(" *********** Begin generating ", numReplicate, " data sets *************\n")
  }
  # for each data set to be generated
  for(b in 1:numReplicate)
  { # get file name
    datafileName<-paste(fileName, "_", b, ".dat", sep="")
    if(!quiet)
    { cat(" *********** data set >> ", datafileName, " *************\n")
      cat(" *** Step 2.1: Generate membership and cluster sizes ***\n") 
    }
    # total number of dimension
    p<-numNonNoisy+numNoisy
    # generate membership and cluster sizes
    tmpmem<-genMemSize(clustszind, numClust, clustSizeEq, 
                       rangeN, clustSizes, p, quiet)
    # mem is the set of the memberships of data points
    mem<-tmpmem$mem
    # size are the numbers of data points in clusters
    size<-tmpmem$size
    # N is the total number of data points
    N<-tmpmem$N 
    if(!quiet)
    { cat(" *********** data set >> ", datafileName, " *************\n")
      cat(" *** Step 2.2: Generate covariance matrices ***\n") 
    }
    # generate covariance matrices
    tmpms<-genMeanCov(numNonNoisy,numClust,
                      rangeVar, sepVal, lambdaLow,
                      ratioLambda, covMethod, alphad, eta, 
                      iniProjDirMethod, projDirMethod, 
                      alpha, ITMAX, eps, quiet) 
    #generate mean and cov
    # 'thetaMat' is a 'numClust' by 'numNonNoisy' matrix 
    # i-th row is the mean vector for the i-th cluster
    thetaMat<-tmpms$thetaMat 
    # an array of covariance matrices; s[,,i] is the covariance matrix
    # for the i-th cluster
    s<-tmpms$s
    # The scalar 'A' is for the length of the simplex, 
    # the vertexes of which allocate
    # cluster centers. The scalar 'A' is to adjust the distance between
    # two clusters so that the separation of the two clusters with 
    # specified covariance matrices is equal to 'sepVal'.
    A<-tmpms$A
    # 'neighbors.mat' is a 'numClust' by 6 matrix. Each row corresponds to a cluster.
    # Columns
    #  1 -- cluster label
    #  2 -- cluster label of its nearest neighboring cluster
    #  3 -- separation index for the cluster and its nearest neighboring cluster
    #  4 -- cluster label of its farthest neighboring cluster
    #  5 -- separation index for the cluster and its farthest neighboring cluster
    #  6 -- median separation index for the cluster and its neighboring clusters
    neighbors.mat<-tmpms$neighbors.mat
    # 'egvaluesMat' is a 'numClust' by 'numNonNoisy' matrix. Each row correspond a cluster.
    # rows are the eigenvalues of covariance matrices of the cluster
    egvaluesMat<-tmpms$egvaluesMat
    # cluster fractions
    mypi<-size/sum(size, na.rm=TRUE)
    # Obtain the mean vector and covariance matrix for noisy variables
    tmp<-genNoisyMeanCov(numNoisy, mypi, numClust, numNonNoisy, 
                         thetaMat, s, covMethod, alphad, eta, rangeVar, 
                         eps, quiet)
    # 'thetaMat' is a 'numClust' by 'p' matrix 
    # i-th row is the mean vector for the i-th cluster
    thetaMat<-tmp$thetaMat 
    # an array of covariance matrices; s[,,i] is the covariance matrix
    # for the i-th cluster
    s<-tmp$s  
    # number of variables = numNonNoisy+numNoisy 
    p<-tmp$p 
    # update mean vectors and covariance matrices of variables
    # e.g. rotating data and randomizing the order of variables
    tmp<-updateMeanCov(rotateind, thetaMat, s, numClust, p, 
                       numNoisy, quiet)
    # 'muMat' is a 'numClust' by 'p' matrix 
    # i-th row is the mean vector for the i-th cluster
    muMat<-tmp$muMat
    # an array of covariance matrices; SigmaArray[,,i] is the covariance matrix
    # for the i-th cluster
    SigmaArray<-tmp$SigmaArray
    # set of noisy variables
    noisySet<-tmp$noisySet
    # Q is the rotation matrix (orthogonal matrix)
    Q<-tmp$Q; 
    if(!quiet)
    { cat(" *********** data set >> ", datafileName, " *************\n")
      cat(" *** Step 2.3: Generate data set ***\n") 
    }
    # generate data based on the mean vectors, covariance matrices
    # and cluster sizes.
    y<-genData(mem, N, numClust, p, muMat, SigmaArray, size, quiet)

    if(!quiet)
    { cat(" *********** data set >> ", datafileName, " *************\n")
      cat(" *** Step 2.4: Get separation indices and projection directions ***\n") 
    }
    # get the theoretical separation index and projection direction
    tmpsep<-getSepProjTheory(muMat, SigmaArray, iniProjDirMethod, 
                             projDirMethod, alpha, ITMAX, eps, quiet)
    # 'sepValMat' is a 'numClust' by 'numClust' matrix
    # 'sepValMat[i,j]' is the separation index between clusters i and j
    sepValMat<-tmpsep$sepValMat
    # 'myprojDir' is a 'numClust' by 'numClust' by 'p' array
    # 'myprojDir[i,j,]' is the projection direction for clusters i and j
    myprojDir<-tmpsep$projDirArray
    if(outputEmpirical)
    { # get empirical separation indices and projection directions
      tmp<-getSepProjData(y, mem, iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet)
      # 'Jhat2' is a 'numClust' by 'numClust' matrix
      # 'Jhat2[i,j]' is the separation index between clusters i and j
      Jhat2<-tmp$sepValMat
      # 'empProjDir' is a 'numClust' by 'numClust' by 'p' array
      # 'empProjDir[i,j,]' is the projection direction for clusters i and j
      empProjDir<-tmp$projDirArray
    } 
    else { 
      empProjDir<-NULL 
      Jhat2<-NULL
    }

    if(!quiet)
    { cat(" *********** data set >> ", datafileName, " *************\n")
      cat(" *** Step 2.5: Generating outliers ***\n") 
    }
    # Generate outliers
    tmpout<-genOutliers(numOutlier, y)
    # number of outliers
    nOut<-tmpout$nOut
    # 'y.out' is a 'nOut' by 'p' matrix, where 'nOut' is the number of outliers
    y.out<-tmpout$y.out
    if(nOut>0)
    { # add outliers
      y<-rbind(y, y.out)
      # the memberships of outliers are zero
      mem<-c(mem, rep(0, nOut))
    }

    if(!quiet)
    { cat(" *********** data set >> ", datafileName, " *************\n")
      cat(" *** Step 2.6: Output log information ***\n") 
    }
    # output log information, e.g. mean vectors and covariance matrices, etc.
    if(outputLogFlag)
    { outputLog(b, fileName, alpha, sepVal, numClust, size,
        N, p, numNoisy, noisySet, nOut, thetaMat, s, muMat, SigmaArray, Q, 
        sepValMat, Jhat2, myprojDir, empProjDir, egvaluesMat, 
        quiet, outputEmpirical) 
    }
    # output data set

    # output membership
    if(outputDatFlag)
    { memfileName<-paste(fileName, "_", b, ".mem", sep="") #membership
      write.table(mem,file=memfileName,quote=FALSE,
                row.names=FALSE,col.names=FALSE)
    }

    # output noisy variables
    if(outputDatFlag)
    { noisyfileName<-paste(fileName, "_",b,".noisy",sep="") #noisy var
      write.table(t(noisySet),file=noisyfileName,
                  quote=FALSE,row.names=FALSE,
                  col.names=FALSE)
    }

    # record theoretical and sample separation indices
    infoMatTheory<-neighbors.mat
    if(outputEmpirical)
    { # for empirical values, we get the corresponding 'numClust' by '6' matrix
      infoMatData<-nearestNeighborSepVal(Jhat2) 
    } else {
      infoMatData<-NULL
    }

    nn<-nrow(infoMatTheory)
    
    if(b>1)  
    { infoFrameTheory<-rbind(infoFrameTheory, infoMatTheory)
      if(outputEmpirical)
      { infoFrameData<-rbind(infoFrameData, infoMatData) }
      else { infoFrameData<-NULL }
      fileNameVec<-c(fileNameVec, rep(paste(fileName, "_", b, sep=""), nn))

      datList[[b]]<-outputData(b, fileName, y, p, outputDatFlag) 
      memList[[b]]<-mem 
      noisyList[[b]]<-noisySet 
    } else { # the first data set
      infoFrameTheory<-infoMatTheory
      if(outputEmpirical)
      { infoFrameData<- infoMatData }
      else { infoFrameData<-NULL }
      fileNameVec<-rep(paste(fileName, "_", b, sep=""), nn)
      datList<-list(outputData(b, fileName, y, p, outputDatFlag)) 
      memList<-list(mem)
      noisyList<-list(noisySet)
    } 
 
  }
  infoFrameTheory<-data.frame(cluster=infoFrameTheory[,1],
                      nearestClust=infoFrameTheory[,2],
                      nearestSep=infoFrameTheory[,3], 
                      farthestClust=infoFrameTheory[,4],
                      farthestSep=infoFrameTheory[,5], 
                      medianSep=infoFrameTheory[,6])
  infoFrameTheory$fileName<-fileNameVec
  if(outputEmpirical)
  { infoFrameData<-data.frame(cluster=infoFrameData[,1],
                      nearestClust=infoFrameData[,2],
                      nearestSep=infoFrameData[,3], 
                      farthestClust=infoFrameData[,4],
                      farthestSep=infoFrameData[,5], 
                      medianSep=infoFrameData[,6])
    infoFrameData$fileName<-fileNameVec
  } else {
    infoFrameData<-NULL
  }

  tmpnames<-paste(fileName, "_", 1:numReplicate, sep="")
  names(datList)<-tmpnames
  names(memList)<-tmpnames
  names(noisyList)<-tmpnames

  if(outputInfo)
  { fileNameInfo<-paste(fileName, "_info.log", sep="")
    msg<-"Theoretical separation information data frame>>>>>>>>>>>\n"
    write.table(msg, append=FALSE, file=fileNameInfo, quote=FALSE,
                row.names=FALSE, col.names=FALSE)
    msg<-colnames(infoFrameTheory)
    write(msg, append=TRUE, file=fileNameInfo, ncolumns=7)
    tt<-infoFrameTheory
    tt[,1:6]<-round(tt[,1:6], 3)
    write.table(tt, file=fileNameInfo, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE, append=TRUE)

    msg<-"\nEmpirical separation information data frame>>>>>>>>>>>\n"
    write.table(msg, append=TRUE, file=fileNameInfo, quote=FALSE,
                row.names=FALSE, col.names=FALSE)
    msg<-colnames(infoFrameData)
    write(msg, append=TRUE, file=fileNameInfo, ncolumns=7)
    if(outputEmpirical)
    { tt<-infoFrameData
      tt[,1:6]<-round(tt[,1:6], 3)
      write.table(tt, file=fileNameInfo, quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE, append=TRUE)
    }
  }
  if(!quiet)
  { 
    cat(" *********** End of generating ", numReplicate, " data sets *************\n")
  }

  invisible(list(infoFrameTheory=infoFrameTheory, infoFrameData=infoFrameData,
                 datList=datList, memList=memList, noisyList=noisyList))
}

# Get numNonNoisy+1 vertexs of a simplex in numNonNoisy dimension
# numNonNoisy is the non-noisy number of dimension
genVertexes<-function(numNonNoisy)
{ numNonNoisy<-as.integer(numNonNoisy) 
  if(numNonNoisy<=1 || !is.integer(numNonNoisy))
  {
    stop("The number of non-noisy variables should be a positive integer that is greater than 1!\n")
  }
  mx<-numNonNoisy+1; # the number of vertices
  # The i-th row of the matrix "vert" is the coordinates of the i-th vertex
  vert<-matrix(0,nrow=mx, ncol=numNonNoisy)
  vert[1,1]<- -1.
  vert[2,1]<-1.
  for(p in 3:mx)
  { # Value of vert[p] to keep constant distance (of 2) between vertices 
    # get centre of previous vertices, then add co-ordinate in pth dimension
    dd<-0.0
    for(j in 1:(p-2))
    { s<-sum(vert[1:(p-1),j], na.rm=TRUE)
      tem<-s/(p-1.)
      vert[p,j]<-tem
      if(j==1) tem<-1. 
      dd<-dd+tem*tem
    }
    vert[p,p-1]<-sqrt(4.-dd)
  }
  return(vert)
}

# Get G vertices in numNonNoisy dimension, where G can be any integer > 0
genShiftedVertexes<-function(G, numNonNoisy)
{ 
  G<-as.integer(G)
  if(G<1 || !is.integer(G))
  {
    stop("The number of vertices 'G' should be a positive integer!\n")
  }
  numNonNoisy<-as.integer(numNonNoisy)
  if(numNonNoisy<=1 || !is.integer(numNonNoisy))
  {
    stop("The number of non-noisy variables should be a positive integer that is greater than 1!\n")
  }
 
  # First get numNonNoisy+1 vertices in a numNonNoisy-dimensional space
  vertex<-genVertexes(numNonNoisy)
  if(G<=numNonNoisy+1) 
  { 
    return(vertex[1:G,]) 
  }
  #For G>numNonNoisy+1, start adding vertices from the following sequence after v(numNonNoisy+1).
  #v(2)+2*e1, ..., v(numNonNoisy+1)+2*e1,
  #v(2)+4*e1, ..., v(numNonNoisy+1)+4*e1,
  #v(2)+6*e1, ..., v(numNonNoisy+1)+6*e1,
  #Essentially this just keeps on adding points on a shifted symmetric simplex.
  # e1 is the numNonNoisy x 1 vector whose elements are all zero except that e1[1]=1.
  vertex2<-matrix(0, nrow=G, ncol=numNonNoisy)
  vertex2[1:(numNonNoisy+1),]<-vertex
  e1<-rep(0,numNonNoisy)
  e1[1]<-1
  tmpNumNonNoisy<-(G-numNonNoisy-1)/numNonNoisy
  nG<-floor(tmpNumNonNoisy)
  res<-(G-numNonNoisy-1)%%numNonNoisy 
  m<-numNonNoisy+1
  # m indicates the label of the current cluster
  # nG indicates how many shifted symmetric simplex we need
  if(nG>0)
  { for(j in 1:nG) 
    { for(i in 1:numNonNoisy) 
      { m<-m+1
        vertex2[m,]<-vertex[i+1,]+2*j*e1
      }
    }
  }
  # res is the number of remindar vertices.
  if(res>0)
  { for(i in 1:res) 
    { m<-m+1
      vertex2[m,]<-vertex[i+1,]+2*(nG+1)*e1 
    }
  }
  rownames(vertex2)<-paste("cluster", 1:G, sep="")
  colnames(vertex2)<-paste("variable", 1:numNonNoisy, sep="")
  return(vertex2) 
}

# Generate mean vectors and covariance matrices for non-noisy variables
# numNonNoisy -- number of non-noisy variables
# G -- number of clusters
# See documentation of genRandomClust for explanation of arguments:
# rangeVar, sepVal, lambdaLow, ratioLambda, covMethod, iniProjDirMethod, 
# projDirMethod, alpha, ITMAX, eps, quiet
genMeanCov<-function(numNonNoisy, G, rangeVar, sepVal, lambdaLow=1, 
                     ratioLambda=10, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
                     alphad=1, eta=1, iniProjDirMethod=c("SL", "naive"), 
                     projDirMethod=c("newton", "fixedpoint"),
                     alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  covMethod<-match.arg(arg=covMethod, choices=c("eigen", "onion", "c-vine", "unifcorrmat"))
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  numNonNoisy<-as.integer(numNonNoisy)

  # check for valid inputs
  if(numNonNoisy<2 || !is.integer(numNonNoisy))
  {
    stop("The number 'numNonNoisy' of non-noisy variables should be a positive integer greater than 1!\n")
  }
  G<-as.integer(G)
  if(G<1 || !is.integer(G))
  {
    stop("The number 'G' of clusters should be positive integer!\n")
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
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  if(sepVal<= -0.999 || sepVal >= 0.999) 
  { 
    stop("The desired separation index 'sepVal' should be in the range (-1, 1)!\n")
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
  # end of checks of valid inputs

  s<-array(0, c(numNonNoisy,numNonNoisy,G))
  dimnames(s)<-list(NULL, NULL, paste("cluster", 1:G, sep=""))
  egvaluesMat<-matrix(0, nrow=G, ncol=numNonNoisy)
  rownames(egvaluesMat)<-paste("cluster", 1:G, sep="")
  colnames(egvaluesMat)<-paste("variable", 1:numNonNoisy, sep="")
  if(G==1) # only one cluster
  { thetaMat<-matrix(0, nrow=1, ncol=numNonNoisy)
    rownames(thetaMat)<-"cluster1"
    colnames(thetaMat)<-paste("variable", 1:numNonNoisy, "\n")
    tmp<-genPositiveDefMat(numNonNoisy, covMethod, alphad, eta, rangeVar, 
         lambdaLow, ratioLambda)
    s[,,1]<-tmp$Sigma
    egvaluesMat[1,]<-tmp$egvalue
    # The scalar 'A' is for the length of the simplex, 
    # the vertexes of which allocate
    # cluster centers. The scalar 'A' is to adjust the distance between
    # two clusters so that the separation of the two clusters with 
    # specified covariance matrices is equal to 'sepVal'.
    A<-0
    neighbors.mat<-NULL
    if(!quiet)
    { cat("Warning: only one cluster in function 'genMeanCov'!\n 'neighbors.mat' is set to be 'NULL'\n")
    }
  } 
  else 
  { a<-rep(0,numNonNoisy)
    a[1]<-1
    asa<-rep(0, G); # sqrt(a^T s a)
    for(i in 1:G) # generate covariance matrices
    { tmp<-genPositiveDefMat(numNonNoisy, covMethod, alphad, eta, rangeVar, lambdaLow, ratioLambda) 
      s[,,i]<-tmp$Sigma
      asa[i]<-sqrt(as.vector(t(a)%*%s[,,i]%*%a)) # the sd of projected data
      egvaluesMat[i,]<-tmp$egvalues
    }
    # obtain the largest two standard deviations
    tmpsd<-sort(asa, decreasing=TRUE)
    s1<-sqrt(tmpsd[1])
    s2<-sqrt(tmpsd[2])
    if(!quiet) 
    { 
      cat(" *** Step 2.2:  Get the suitable value of A.***\n")
    }
    # calculate the initial upper bound of A
    # The scalar 'A' is for the length of the edge of the simplex, 
    # the vertexes of which allocate cluster centers. The scalar 
    # 'A' is to adjust the distance between
    # two clusters so that the separation of the two clusters with 
    # specified covariance matrices is equal to 'sepVal'.
    za<-qnorm(1-alpha/2)
    A2<-(1+sepVal)*za*(s1+s2)/(2*(1-sepVal))
    A2<-A2+1
    for(gg in 1:G)
    { eg<-eigen(s[,,gg])$values
    }
    tmp<-getA2(G, s, sepVal, A2, 
        iniProjDirMethod, projDirMethod, alpha, ITMAX, eps) # get value of A
    if(!quiet)
    { cat(" *** Step 2.3:  Generate mean vectors for clusters ***\n") }
    A<-tmp$minA
    thetaMat<-A*tmp$vertex # obtain the centers of the clusters

    tmp<-getSepProjTheory(thetaMat, s, iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet)
    sepValMat<-tmp$sepValMat
    #projDirArray<-tmp$projDirArray
    d<-as.matrix(dist(thetaMat))
    neighbors.mat<-nearestNeighbor(d, A, sepValMat)
    if(!quiet)
    { # output intermediate results
      cat("true sepVal=", sepVal, "\n")
      cat("before scaling the covariance matrices, neighbors.mat>>\n");
      cat("1st column has the labels of clusters\n");
      cat("2nd column has the labels of its nearest neighboring cluster\n");
      cat("3rd column has the separation indices of the clusters to their nearest neighbors\n")
      cat("4nd column has the labels of its farthest neighboring cluster\n")
      cat("5rd column has the separation indices of the clusters to their farthest neighbors\n")
      cat("6th column has the median separation indices of the clusters to their neighbors\n")
      print(neighbors.mat);
    }

    # Refine Covariance matrices so that the separation indices between
    # any cluster and its nearest neighboring cluster is equal to sepVal.
    tmp<-refineCov(thetaMat, s, A, sepVal, iniProjDirMethod, projDirMethod, 
                   alpha, ITMAX, eps, quiet)
    s<-tmp$s

    tmp<-getSepProjTheory(thetaMat, s, iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet)
    sepValMat<-tmp$sepValMat
    d<-as.matrix(dist(thetaMat))
    neighbors.mat<-nearestNeighbor(d, A, sepValMat)
    if(!quiet)
    { # output intermediate results
      # recalcuate the separation index matrix
      cat("true sepVal=", sepVal, "\n")
      cat("before scaling the covariance matrices, neighbors.mat>>\n");
      cat("1st column has the labels of clusters\n");
      cat("2nd column has the labels of its nearest neighboring cluster\n");
      cat("3rd column has the separation indices of the clusters to their nearest neighbors\n")
      cat("4nd column has the labels of its farthest neighboring cluster\n")
      cat("5rd column has the separation indices of the clusters to their farthest neighbors\n")
      cat("6th column has the median separation indices of the clusters to their neighbors\n")
      print(neighbors.mat);
    }
  }
  return(list(thetaMat=thetaMat, s=s, A=A, neighbors.mat=neighbors.mat,
              egvaluesMat=egvaluesMat))
}

# Refine Covariance matrices so that the separation indices between
# clusters and their nearest neighboring clusters are equal to sepVal.
# thetaMat -- cluster center matrix, 
#    thetaMat[i,] is the cluster center for the i-th cluster.
# s -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# A -- scalar so that the minimum separation index among clusters is equal
#      to sepVal.
#      The scalar 'A' is for the length of the simplex, 
#      the vertexes of which allocate cluster centers. 
#      The scalar 'A' is used to adjust the distance between
#      two clusters so that the separation of the two clusters with 
#      specified covariance matrices is equal to 'sepVal'.
# sepVal -- the minimum separation index set as a priori
# See documentation of genRandomClust for explanation of arguments:
# iniProjDirMethod, projDirMethod, alpha, ITMAX, eps
refineCov<-function(thetaMat, s, A, sepVal, 
                    iniProjDirMethod=c("SL", "naive"), 
                    projDirMethod=c("newton", "fixedpoint"), 
                    alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(A<0)
  { 
    stop("The scalar 'A' should be positive!\n")
  }
  if(sepVal<= -0.999 || sepVal >= 0.999) 
  { 
    stop("The desired separation index 'sepVal' should be in the range (-1, 1)!\n")
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

  k.max.old<-0
  myc<-1
  while(1) 
  { # We first find neighboring clusters for each cluster. Then we calculate
    # the separation index of each cluster to its nearest neighboring
    # cluster. Then we choose the cluster whose minimum separation index is
    # the largest.
    tmp<-neighborsMaxMin(thetaMat, s, A, iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet)
    k.maxmin<-tmp$k.maxmin
    neighbors.maxmin<-tmp$neighbors.maxmin 
    if(k.maxmin==k.max.old) { break }
    k.max.old<-k.maxmin
    maxminJ<-tmp$maxminJ
    # If maxminJ is already equal to sepVal, then we do not need rescale
    # covariance matrices further.
    if(abs(maxminJ-sepVal)<eps) { break }

    # rescale the covariance matrix of the cluster "k.maxmin" so that the
    # separation index between the cluster "k.maxmin" and its nearest
    # neighboring clusters is as close as to sepVal as possible.
    pos<-c(k.maxmin, neighbors.maxmin)
    theta<-thetaMat[pos,]
    Sigma<-s[,,pos]

    # search the upper bound of the myc so that the minimum separation index
    # after scaling is smaller than the true value sepVal
    myupper<-findUpperC(theta, Sigma, sepVal, iniProjDirMethod,
                        projDirMethod, alpha, ITMAX, eps, quiet)
    # find a suitable scalar myc such that the minimum separation index of
    # the cluster k.maxmin is equal to sepVal.
    newfit<-0
    class(newfit)<-"try-error"
    newfit<-try(
         tmp<-uniroot(f=diffSepVal, lower=0.9, upper=myupper, 
      thetaMat=theta, s=Sigma, sepVal=sepVal, 
      iniProjDirMethod=iniProjDirMethod, projDirMethod=projDirMethod, 
      alpha=alpha, ITMAX=ITMAX, eps=eps, quiet=quiet)
    )
    if(sum(class(newfit)=="try-error", na.rm=TRUE))
    { warning("Could not find suitable upper bound of 'myc'!\n 'myc' is set to be 'myupper'!\n")
      myc <-myupper 
    } else {
      myc<-tmp$root
    }
    # scale the covariance matrix of the cluster k.maxmin so that the
    # minimum separation index of the cluster k.maxmin is equal to sepVal.
    s[,,k.maxmin]<-myc*s[,,k.maxmin]
  }
  return(list(s=s, myc=myc))
}

# We first find neighboring clusters for each cluster. Then we calculate
# the separation index of each cluster to its nearest neighboring
# cluster. Then we choose the cluster whose minimum separation index is
# the largest.
# thetaMat -- cluster center matrix. thetaMat[i,] is the cluster center
#             for the i-th cluster
# s -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# A -- scalar so that the minimum separation index among clusters is equal
#      to sepVal
#     The scalar 'A' is for the length of the simplex, 
#     the vertexes of which allocate cluster centers. 
#     The scalar 'A' is used to adjust the distance between
#     two clusters so that the separation of the two clusters with 
#     specified covariance matrices is equal to 'sepVal'.
# See documentation of genRandomClust for explanation of arguments:
# iniProjDirMethod, projDirMethod, alpha, eps, ITMAX, quiet
neighborsMaxMin<-function(thetaMat, s, A, 
                          iniProjDirMethod=c("SL", "naive"), 
                          projDirMethod=c("newton", "fixedpoint"), 
                          alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(A<0)
  { 
    stop("The scalar 'A' should be positive!\n")
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

  d<-as.matrix(dist(thetaMat))
  numClust<-nrow(thetaMat)
  p<-ncol(thetaMat)
  # calculate the separation index matrix
  tmp<-getSepProjTheory(thetaMat, s, iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet)
  sepValMat<-tmp$sepValMat
  # Find neighboring clusters of each cluster and find the cluster whose
  # minimum separation index with its neighboring clusters is the largest.
  tmp<-findNeighbors(d, A, sepValMat)
  k.maxmin<-tmp$k.maxmin; # the cluster which we are interested in
  maxminJ<-tmp$maxminJ; # the maximum minimum separation index
  neighbors<-tmp$neighbors # neighboring clusters of each cluster
  # find the neighboring clusters of the cluster k.maxmin, whose separation
  # indices with the cluster k.maxmin is equal to maxminJ
  tmpJ<-sepValMat[k.maxmin,]
  k2<-which(tmpJ==maxminJ)
  # To check if the neighboring clusters in the set k2 have larger minimum
  # separation indices with their neighboring clusters than maxminJ
  len<-length(k2)
  minJ<-rep(0,len)
  for(i in 1:len)
  { pos<-k2[i]
    tmpJ<-sepValMat[pos,]
    tmpJ[pos]<-1
    minJ[i]<-min(tmpJ, na.rm=TRUE)
  }
  maxminJ2<-max(minJ, na.rm=TRUE)
  pos<-which(minJ==maxminJ2)
  pos<-k2[pos[1]]
  if(maxminJ<maxminJ2)
  { k.maxmin<-pos; maxminJ<-maxminJ2; }
  neighbors.maxmin<-neighbors[[k.maxmin]]
  res<-list(k.maxmin=k.maxmin, maxminJ=maxminJ, neighbors.maxmin=neighbors.maxmin) 
  return(res)
}

# Find a crude upper bound 'myc' such that after scaling the covariance 
# matrix of the first cluster, the minimum separation index of first cluster
# is less than sepVal. Before scaling the minimum separation index is larger
# than sepVal. This is why we need to scale the covariance matrix of the
# first cluster.
# thetaMat -- cluster center matrix. 
# s -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# sepVal -- the minimum separation index set as a priori
# See documentation of genRandomClust for explanation of arguments:
# alpha, iniProjDirMethod, projDirMethod, ITMAX, eps, quiet
findUpperC<-function(thetaMat, s, sepVal, 
                     iniProjDirMethod=c("SL", "naive"), 
                     projDirMethod=c("newton", "fixedpoint"), 
                     alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  if(sepVal<= -0.999 || sepVal >= 0.999) 
  { 
    stop("The desired separation index 'sepVal' should be in the range (-0.999, 0.999)!\n")
  }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<=0 || !is.integer(ITMAX))
  {
    stop("The maximum iteration number allowed 'ITMAX' should be a positive integer!\n")
  }

  myc<-2
  numClust<-nrow(thetaMat) 
  p<-ncol(thetaMat)
  while(1)
  { tmps<-s
    # scale the covariance matrix of the first cluster
    tmps[,,1]<-myc*s[,,1]
    # calculate the separation index matrix
    tmp<-getSepProjTheory(thetaMat, tmps, iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet)
    sepValMat<-tmp$sepValMat
    sepVal1<-sepValMat[1,]
    sepVal1[1]<-1
    # find the minimum separation index of the first cluster
    minsepVal<-min(sepVal1, na.rm=TRUE)
    if(minsepVal<sepVal) { break }
    myc<-myc^2
  }
  return(myc)
}

# Calculate the difference between the minimum separation index and sepVal.
# The minimum separation index is between the cluster 1 and its
# nearest neighbor. The covariance matrix of the cluster 1 is scaled by a
# multiplier "myc".
# This function is used as the objective function to obtain the optimal
# scalar myc.
# myc -- the scalar
# thetaMat -- cluster center matrix. thetaMat[i,] is the cluster center
#      for the i-th cluster
# s -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# sepVal -- the minimum separation index set as a priori
# See documentation of genRandomClust for explanation of arguments:
# alpha, iniProjDirMethod, projDirMethod, ITMAX, eps, quiet
diffSepVal<-function(myc, thetaMat, s, sepVal, 
                    iniProjDirMethod=c("SL", "naive"), 
                    projDirMethod=c("newton", "fixedpoint"), 
                    alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(myc<0)
  {
    stop("The scalar 'myc' should be positive!\n")
  }
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  if(sepVal<= -0.999 || sepVal >= 0.999) 
  { 
    stop("The desired separation index 'sepVal' should be in the range (-0.999, 0.999)!\n")
  }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<=0 || !is.integer(ITMAX))
  {
    stop("The maximum iteration number allowed 'ITMAX' should be a positive integer!\n")
  }

  tmps<-s
  tmps[,,1]<-myc*tmps[,,1]
  numClust<-nrow(thetaMat)
  p<-ncol(thetaMat)
  tmp<-getSepProjTheory(thetaMat, tmps, iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet)
  sepValMat<-tmp$sepValMat
  sepVal1<-sepValMat[1,]
  sepVal1[1]<-1
  minsepVal<-min(sepVal1, na.rm=TRUE)
  res<-minsepVal-sepVal
  return(res)
}

# Find the neighboring cluster of the cluster k.
# dk is set of the distances from cluster k to other clusters
# dk[k,k]=0, i.e. the distance from cluster k to itself is zero.
# A -- scalar so that the minimum separation index among clusters is equal
#      to sepVal
myNeighbor<-function(d, k, A)
{ 
  k<-as.integer(k)
  if(k<0 || !is.integer(k))
  {
    stop("The cluster label 'k' should be positive integer!\n")
  }
  if(A<0)
  {
    stop("The scalar 'A' should be positive\n")
  }
  len<-length(d)
  if(len<1)
  { 
    stop("The distance vector 'd' should contains at least 1 element!\n")
  }
  pos<-which(d<0)
  if(length(pos)>0)
  {
      stop("The element of the distance vector 'd' should be positive!\n")
  }

  dk<-as.vector(d[k,])
  # In our design, the distances from a cluster to its neighboring clusters
  # are equal to A*2.
  tmp<-which(abs(dk-A*2)<1.0e-10)
  tmp<-sort(tmp)
  return(tmp)
}

# Find neighboring clusters of each cluster and find the cluster whose
# minimum separation index with its neighboring clusters is the largest.
# d -- distance matrix of numClust clusters
# A -- scalar so that the minimum separation index among clusters is equal
#      to sepVal
# sepValMat -- the separation index matrix
findNeighbors<-function(d, A, sepValMat)
{ 
  tmp<-as.vector(d)
  len<-length(tmp)
  if(len<1)
  { 
    stop("The distance matrix 'd' should contains at least 1 element!\n")
  }
  pos<-which(tmp<0)
  if(length(pos)>0)
  {
      stop("The element of the distance matrix 'd' should be positive!\n")
  }
  if(A<0)
  {
    stop("The scalar 'A' should be positive\n")
  }
  tmp<-as.vector(sepValMat)
  len<-length(tmp)
  if(len<1)
  { 
    stop("The separation index matrix 'sepValMat' should contains at least 1 element!\n")
  }
  pos<-which(tmp< -1 || tmp>1)
  if(length(pos)>0)
  {
      stop("The element of the separation index matrix 'sepValMat' should be between [-1, 1]!\n")
  }

  numClust<-nrow(sepValMat)
  # find the neighboring clusters of the first cluster
  neighbors<-list(myNeighbor(d, 1, A))
  # the separation indices between 1st cluster and its neighboring clusters
  sepVal.neighbors<-list(sepValMat[1,neighbors[[1]]])
  tmpJ<-as.vector(sepVal.neighbors[[1]])
  maxminJ<-min(tmpJ, na.rm=TRUE)
  k.maxmin<-1
  for(k in 2:numClust)
  { neighbors[[k]]<-myNeighbor(d, k, A) 
    sepVal.neighbors[[k]]<-sepValMat[k,neighbors[[k]]]
    tmpJ<-as.vector(sepVal.neighbors[[k]])
    tmpminJ<-min(tmpJ, na.rm=TRUE)
    if(tmpminJ>maxminJ) 
    { k.maxmin<-k 
      maxminJ<-tmpminJ
    } 
  }
  res<-list(k.maxmin=k.maxmin, maxminJ=maxminJ, neighbors=neighbors, 
            sepVal.neighbors=sepVal.neighbors) 
  return(res)
}

# Separation information matrix containing
# the nearest neighbor and farthest neighbor of each cluster
# sepValMat -- the separation index matrix
nearestNeighborSepVal<-function(sepValMat)
{ 
  if(!is.matrix(sepValMat))
  {
    stop("'sepValMat' should be a matrix!\n")
  }
  if(nrow(sepValMat) != ncol(sepValMat))
  {
    stop("'sepValMat' should be a squared matrix!\n")
  }

  tmp<-as.vector(sepValMat)
  len<-length(tmp)
  if(len<1)
  { 
    stop("The separation index matrix 'sepValMat' should contains at least 1 element!\n")
  }
  pos<-which(tmp< -1 || tmp>1)
  if(length(pos)>0)
  {
      stop("The element of the separation index matrix 'sepValMat' should be between [-1, 1]!\n")
  } 

  numClust<-nrow(sepValMat)

  # 'neighbors.mat' is a 'numClust' by '6' matrix.
  # 1st column has the labels of clusters
  # 2nd column has the labels of its nearest neighboring cluster
  # 3rd column has the separation indices of the clusters to their nearest
  #   neighbors
  # 4nd column has the labels of its farthest neighboring cluster
  # 5rd column has the separation indices of the clusters to their farthest
  #   neighbors
  # 6th column has the median separation indices of the clusters to their
  #   neighbors
  neighbors.mat<-matrix(0, nrow=numClust, ncol=6)
  rownames(neighbors.mat)<-1:numClust
  colnames(neighbors.mat)<-c("cluster", "neighbor_nearest", "sep_nearest", "neighbor_farthest", "sep_farthest", "sep_median")
  neighbors.mat[,1]<-1:numClust
  for(i in 1:numClust)
  { Jvec<-as.vector(sepValMat[i,])
    pos.max<-which(Jvec==max(Jvec))[1]
    neighbors.mat[i,4]<-pos.max
    neighbors.mat[i,5]<-Jvec[pos.max]
    Jvec2<-Jvec[-i]
    minJ<-min(Jvec2, na.rm=TRUE)
    neighbors.mat[i,3]<-minJ
    pos.min<-which(Jvec==minJ)[1]
    neighbors.mat[i,2]<-pos.min
    medianJ<-median(Jvec2, na.rm=TRUE)
    neighbors.mat[i,6]<-medianJ
  }
  return(neighbors.mat)
}

# find the nearest neighbor and farthest neighbor of each cluster
# d -- distance matrix of numClust clusters
# A -- scalar so that the minimum separation index among clusters is equal
#      to sepVal
# sepValMat -- the separation index matrix
nearestNeighbor<-function(d, A, sepValMat)
{ # 1st column has the labels of clusters
  # 2nd column has the labels of its nearest neighboring cluster
  # 3rd column has the separation indices of the clusters to their nearest
  #   neighbors
  # 4nd column has the labels of its farthest neighboring cluster
  # 5rd column has the separation indices of the clusters to their farthest
  #   neighbors
  # 6th column has the median separation indices of the clusters to their
  #   neighbors
  numClust<-nrow(sepValMat)
  neighbors.mat<-matrix(0, nrow=numClust, ncol=6)
  rownames(neighbors.mat)<-1:numClust
  colnames(neighbors.mat)<-c("cluster", "neighbor_nearest", "sep_nearest", "neighbor_farthest", "sep_farthest", "sep_median")
  neighbors.mat[,1]<-1:numClust
  # find neighbors of 1st cluster
  for(k in 1:numClust)
  { neighbork<-myNeighbor(d, k, A)
    tmpJ<-as.vector(sepValMat[k,neighbork])
    maxJ<-max(tmpJ, na.rm=TRUE)
    pos.max<-neighbork[which(tmpJ==maxJ)]
    minJ<-min(tmpJ, na.rm=TRUE)
    pos.min<-neighbork[which(tmpJ==minJ)]
    medianJ<-median(tmpJ, na.rm=TRUE)
    neighbors.mat[k,2]<-pos.min[1]
    neighbors.mat[k,3]<-minJ
    neighbors.mat[k,4]<-pos.max[1]
    neighbors.mat[k,5]<-maxJ
    neighbors.mat[k,6]<-medianJ
  }
  return(neighbors.mat)
}

# Generate membership and cluster sizes
# clustszind -- cluster size indicator. 
#      clustszind=1 indicates that all cluster have equal size.  
#         The size is specified by clustSizeEq.
#      clustszind=2 indicates that the cluster sizes are randomly 
#         generated from the range rangeN.
#      clustszind=3 indicates that the cluster sizes are specified
#         via the vector clustSizes.
#      The default value is 2 so that the generated clusters are more realistic.
# G -- the number of clusters
# clustSizeEq -- if clustszind=1, then this is the constant cluster size 
#      The default value 100 is a reasonable cluster size.
# rangeN -- if clustszind=2, then the cluster sizes are randomly generaged
#         from the range 'rangeN'. The default range is [50, 200] which
#         can produce reasonable variability of cluster sizes.
# clustSizes -- if clustszind=3, then the cluster sizes are specified via
#         clustSizes. An input is required with
#         the default value of 'clustSizes' set as 'NULL'.
# p -- number of variables (non-noisy and noisy variables)
# quiet -- a flag to switch on/off the outputs of intermediate results.
#          The default value is 'TRUE'.
genMemSize<-function(clustszind, G, clustSizeEq, rangeN, clustSizes, p, quiet=TRUE)
{ 
  if(!is.element(clustszind, c(1,2,3))) 
  { 
    stop("Cluster size indicator 'clustszind' should be 1, 2, or 3!\n")
  }
  clustSizeEq<-as.integer(clustSizeEq)
  if(clustSizeEq<2 || !is.integer(clustSizeEq)) 
  { 
    stop("Cluster size 'clustSizeEq' should be an integer greater than 1!\n") 
  }
  if(length(rangeN)!=2) 
  { 
    stop("The range 'rangeN' for cluster size should be a numeric vector of length 2!\n") 
  }
  if(rangeN[1]>rangeN[2]) 
  { 
    stop("First element of 'rangeN' should be smaller than second!\n")
  }
  rangeN<-as.integer(rangeN)
  if(rangeN[1]<2 || !is.integer(rangeN[1])) 
  { 
    stop("The lower bound 'rangeN[1]' for the range of cluster sizes should be an integer greater than 1!\n") 
  }
  if(!is.integer(rangeN[2])) 
  { 
    stop("The upper bound 'rangeN[2]' for the range of cluster size should be integer!\n") 
  }
  if(clustszind==3)
  {
    len<-length(clustSizes)
    if(len!=G || is.null(clustSizes)) 
    {
      stop("The number of elements in 'clustSizes' should be equal the number 'G' of clusters when the value of 'clustszind' is equal to 3!\n")
    }
    clustSizes<-as.integer(clustSizes)
    for(i in 1:len)
    { if(clustSizes[i]<1 || !is.integer(clustSizes[i]))
      { stop(paste("The cluster size for the ", i, "-th cluster should be an integer greater than 1!\n", sep=""))
      }
    }
  }
  p<-as.integer(p)
  if(p<1 || !is.integer(p))
  {
    stop("The number 'p' of variables should be positive integer!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  if(!quiet)
  { cat(" *** Step 2.4:  Generate cluster sizes and membership of each data point ***\n") }
  if(clustszind==1)
  { # N is the total sample size; mem is the set of the memberships of 
    # data points.
    N<-G*clustSizeEq
    mem<-rep(1:G, rep(clustSizeEq, G)) 
    mem<-sample(mem) # randomize the order of the membership
  } 
  else if(clustszind==2)
  { n.low<-rangeN[1]
    n.upp<-rangeN[2]  
    ratio<-n.upp/n.low
    if(n.low<p)
    { n.low<-max(2*p, 50, na.rm=TRUE)
      n.upp<-ratio*n.low
    }
    n.vec<-sample(n.low:n.upp,G,replace=TRUE) #generate cluster sizes
    N<-sum(n.vec, na.rm=TRUE)
    pihat<-n.vec/N
    mem<-sample(1:G, N, replace=TRUE, prob=pihat);
  } 
  else 
  { # cluster sizes are specified by the user.
    mem <- sample(unlist(lapply(1:G, function(x) rep.int(x, times =
        clustSizes[x]))))
    N <- sum(clustSizes, na.rm=TRUE)
    pihat<-clustSizes/N
  }
  size<-rep(0,G) # get final cluster sizes
  for(i in 1:G) 
  { size[i]<-sum(mem==i, na.rm=TRUE) }
  return(list(mem=mem, size=size, N=N))
}

# Adding noisy variables. Noisy variables have the same distribution across
# clusters. We need to make sure the variations of noisy variables are
# similar to those of non-noisy variables.
# If the variations of noisy varialbes are much smaller than those of non-noisy
# variables, then we implicitly downweighted noisy variables.
#
# numNoisy --- the number of noisy variables. The default values of 'numNoisy'
#     and 'numOutlier' are 0 so that we get 'clean' data sets. 
# mypi -- the proportions of cluster sizes. sum(mypi)=1.
# G -- the number of clusters
# numNonNoisy -- the number of non-noisy variables
# thetaMat -- cluster center matrix. thetaMat[i,] is the cluster center
#      for the i-th cluster
# s -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# See documentation of genRandomClust for explanation of arguments:
# covMethod, rangeVar, eps, quiet
genNoisyMeanCov<-function(numNoisy, mypi, G, numNonNoisy, thetaMat, s, covMethod, alphad, eta, rangeVar, eps=1.0e-10, quiet=TRUE)
{ 
  numNoisy<-as.integer(numNoisy)
  if(numNoisy<0 || !is.integer(numNoisy)) 
  { 
    stop("The number 'numNoisy' of noisy variables should be a non-negative integer!\n") 
  }
  tmp<-which(mypi<0)
  if(length(tmp)>0)
  {
    stop("All elements of 'mypi' should be positive!\n")
  }
  tmp<-which(mypi>1)
  if(length(tmp)>0)
  {
    stop("All elements of 'mypi' should be less than 1!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }
  if(abs(sum(mypi, na.rm=TRUE)-1.0)>eps)
  {
    stop("summation of 'mypi' should be equal to 1!\n")
  } 
  G<-as.integer(G)
  if(G<1 || !is.integer(G))
  {
    stop("The number of clusters should be positive integer!\n")
  }
  numNonNoisy<-as.integer(numNonNoisy)
  if(numNonNoisy<2 || !is.integer(numNonNoisy))
  {
    stop("The number of non-noisy variables should be positive integer greater than 1!\n")
  }
  if(rangeVar[1]<0) 
  { 
    stop("The lower bound 'rangeVar[1]' for the range of variances should be positive!\n") 
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  if(numNoisy==0) { return(list(thetaMat=thetaMat, s=s, p=numNonNoisy)) }
  if(!quiet)
  { cat(" *** Step 2.5:  Generate mean vector and covariance matrix for the noisy variables.***\n") }
  p2<-numNoisy # number of noisy variables
  # Obtain the mean vector and covariance matrix of the mixture of
  # distributions sum_{k=1}^{G}mypi[k]*f_k(X).
  tmp<-meanCovMixture(thetaMat, s, mypi)
  mu.noisy<-tmp$mu.mixture
  Sigma.noisy<-tmp$Sigma.mixture
  # Obtain the range of mu.noisy
  range.mu<-range(mu.noisy, na.rm=TRUE)
  # Generate the mean vector of noisy variables
  mu.noisy<-runif(n=p2, min=range.mu[1], max=range.mu[2])
  # Obtain the range of the eigen values of Sigma.noisy
  egvalues<-eigen(Sigma.noisy, symmetric=TRUE)$values
  range.eg<-range(egvalues, na.rm=TRUE)
  # Generate the covariance matrix of noisy variables
  low<-range.eg[1]; upp<-range.eg[2]; 
  # obtain covariance matrices
  tmp<-genPositiveDefMat(p2, covMethod, alphad, eta, rangeVar, low, upp/low) 
  s.noisy<-tmp$Sigma # obtain covariance matrices
  p<-numNonNoisy+p2 # total number of variables
  # update the mean vectors and covariance matirces
  # The mean of noisy variables are zero
  mu.noisy2<-rep(mu.noisy, G)
  mu.mat<-matrix(mu.noisy2, nrow=p2, ncol=G)
  mu.mat<-t(mu.mat)
  tmpmu.mat<-matrix(0,nrow=G, ncol=p)
  tmpmu.mat[,1:numNonNoisy]<-thetaMat
  tmpmu.mat[,(numNonNoisy+1):p]<-mu.mat
  thetaMat<-tmpmu.mat
  tmps<-array(0, c(p,p,G))
  dimnames(tmps)<-list(NULL, NULL, paste("cluster", 1:G, sep=""))
  for(i in 1:G)
  { tmps[1:numNonNoisy,1:numNonNoisy,i]<-s[,,i]
    tmps[(numNonNoisy+1):p, (numNonNoisy+1):p,i]<-s.noisy
  }
  return(list(thetaMat=thetaMat, s=tmps, p=p)) 
}

# Obtain the mean vector and covariance matrix of the mixture of
# distributions sum_{k=1}^{G}mypi[k]*f_k(X).
# thetaMat -- cluster center matrix. thetaMat[i,] is the cluster center
#      for the i-th cluster
# s -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# mypi -- the proportions of cluster sizes. sum(mypi)=1.
meanCovMixture<-function(thetaMat, s, mypi)
{ 
  tmp<-which(mypi<0)
  if(length(tmp)>0)
  {
    stop("All elements of 'mypi' should be positive!\n")
  }
  tmp<-which(mypi>1)
  if(length(tmp)>0)
  {
    stop("All elements of 'mypi' should be less than 1!\n")
  }

  G<-nrow(thetaMat)
  p<-ncol(thetaMat)
  mu.mixture<-as.vector(mypi[1]*thetaMat[1,])
  Sigma.mixture<-mypi[1]*s[,,1]
  for(k in 2:G)
  { mu.mixture<-mu.mixture+as.vector(mypi[k]*thetaMat[k,]) 
    Sigma.mixture<-Sigma.mixture+mypi[k]*s[,,k]
  }
  term2<-rep(0, p)
  for(k1 in 1:(G-1))
  { muk1<-as.vector(thetaMat[k1,])
    for(k2 in (k1+1):G)
    { muk2<-as.vector(thetaMat[k2,])
      tmpk12<-as.vector(muk1-muk2)
      term2<-term2+mypi[k1]*mypi[k2]*tmpk12%*%t(tmpk12)
    }
  }
  Sigma.mixture<-Sigma.mixture+term2
  return(list(mu.mixture=mu.mixture, Sigma.mixture=Sigma.mixture))
}

# Update mean vectors and covariance matrices of non-noisy variables.
# e.g. rotating data and randomizing the order of variables.
# We should not rotate noisy variables in case that noisy variables are
# changed to non-noisy via rotation. e.g. a bivariate data set, the second
# variable is noisy variable and the first variable is non-noisy variable.
# If we rotate data 90 degree, then the second variable is non-noisy while
# the first variable is noisy.
# We rotate data by the transformation Y=QX. So mu_Y=Q*mu_X,
# Sigma_Y=Q*Sigma_X*Q^T.
#
# rotateind-- if rotateind =TRUE, then rotate data so that we may not detect the
#      full cluster structure from scatter plots, otherwise not rotate.
#      By default, 'rotateind=TRUE' to generate more realistic data sets.
# thetaMat -- cluster center matrix. thetaMat[i,] is the cluster center
#      for the i-th cluster
# s -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# G -- the number of clusters
# p -- the number of variables (both non-noisy and noisy variables)
# numNoisy -- the number of noisy variables
# quiet -- a flag to switch on/off the outputs of intermediate results.
#      The default value is 'TRUE'.
updateMeanCov<-function(rotateind, thetaMat, s, G, p, numNoisy, quiet=TRUE)
{ 
  if(!is.logical(rotateind))
  {
    stop("The value of the rotation indicator 'rotateind' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  G<-as.integer(G)
  if(G<1 || !is.integer(G))
  {
    stop("The number 'G' of clusters should be positive integer!\n")
  }
  p<-as.integer(p)
  if(p<1 || !is.integer(p))
  {
    stop("The number 'p' of variables should be positive integer!\n")
  }
  numNoisy<-as.integer(numNoisy)
  if(numNoisy<0 || !is.integer(numNoisy)) 
  { 
    stop("The number 'numNoisy' of noisy variables should be a non-negative integer!\n") 
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  if(!quiet)
  { cat(" *** Step 2.6:  Generate the numNonNoisy x numNonNoisy orthogonal matrix Q ***\n") }
  numNonNoisy<-p-numNoisy # the number of non-noisy variables
  if(rotateind) { Q<-genOrthogonal(numNonNoisy) } 
  else { Q<-diag(numNonNoisy) }
  if(!quiet)
  { cat(" *** Step 2.7:  Obtain the rotated mean vectors and covariance matrices. The noisy variables will not be rotated. ***\n") 
  }
  muMat<-as.matrix(thetaMat)
  muMat[,1:numNonNoisy]<-muMat[,1:numNonNoisy]%*%t(Q)
  SigmaArray<-s
  for(i in 1:G) 
  { SigmaArray[1:numNonNoisy,1:numNonNoisy,i]<-Q %*% s[1:numNonNoisy,1:numNonNoisy,i]%*%t(Q) }
  if(!quiet)
  { cat(" *** Step 2.8:  Randomize the order of variables (including noisy variables)***\n") 
  }
  # randomize columns
  myset1<-1:p
  myset2<-sample(myset1, replace=FALSE)
  if(numNoisy>0) { mynoisy<-which(myset2>numNonNoisy) } 
  else { mynoisy<-0 }
  muMat<-muMat[,myset2]
  for(i in 1:G)
  { SigmaArray[,,i]<-SigmaArray[myset2, myset2, i] }
  return(list(muMat=muMat, SigmaArray=SigmaArray, noisySet=mynoisy, Q=Q))
}

# Generate Data
# mem -- memberships of data points
# N -- the total number of data points
# G -- the number of clusters
# p -- the total number of variables
# muMat -- cluster center matrix. thetaMat[i,] is the cluster center
#             for the i-th cluster
# SigmaArray -- array of covariance matrices. s[,,i] is the covariance matrix
#      for the i-th cluster
# size -- cluster sizes
# quiet -- indicates if the intermediate results should be output.
genData<-function(mem, N, G, p, muMat, SigmaArray, size, quiet=TRUE)
{ 
  if(length(unique(mem))!=G)
  { stop("The number of clusters obtained from the membership vector 'mem' is not equal to the specified number 'G' of clusters!\n")
  }
  G<-as.integer(G)
  if(G<1 || !is.integer(G))
  {
    stop("The number 'G' of clusters should be positive integer!\n")
  }
  p<-as.integer(p)
  if(p<1 || !is.integer(p))
  {
    stop("The number 'p' of variables should be positive integer!\n")
  }
  if(length(size)!=G)
  {
    stop("The length of 'size' is not equal to the number of clusters!\n")
  }
  tmp<-which(size<1)
  if(length(tmp))
  { stop("cluster sizes should be positive!\n") }
  size<-as.integer(size)
  for(i in 1:length(size))
  { if(!is.integer(size[i]))
    { stop("cluster sizes should be integer!\n") }
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  if(!quiet)
  { cat(" *** Step 2.9:  Generate data set according to the mean vectors and covariance matrices, membership, and cluster sizes.***\n") }
  y<-matrix(0, nrow=N, ncol=p)
  for(i in 1:G)
  { y[mem==i,]<-mvrnorm(size[i], mu=muMat[i,], Sigma=SigmaArray[,,i]) }
  return(y)
}

# Generate outliers. Outliers will be generated from a uniform distribution.
# numOutlier -- number/ratio of outliers. If numOutlier is a positive integer,
#   then numOutlier means the number of outliers. If numOutlier is a real
#   number between (0, 1), then numOutlier means the ratio of outliers, i.e. 
#   the number of outliers is equal to round(numOutlier*N), where N is the 
#   total number of non-outliers. If numOutlier is a real number greater than
#   1, then we round numOutlier to an integer. 
# y -- Nxp data matrix
genOutliers<-function(numOutlier, y)
{ 
  if(numOutlier<0)
  { stop("Number of outliers should be positive!\n") }

  N<-nrow(y)
  p<-ncol(y)
  # obtain the number of outliers
  nOut<-outlierSize(numOutlier, N)
  if(nOut>0)
  { # Find the cooridinate-wise mean and sd of the data matrix y
    mu<-apply(y, 2, mean, na.rm=TRUE)
    s<-apply(y, 2, sd, na.rm=TRUE)
    # Calculate the range of the uniform distribution in each dimension.
    # [mu[j]-5*s[j], mu[j]+5*s[j]], j=1,...,p.
    y.out<-matrix(0,nrow=nOut, ncol=p)
    for(j in 1:p)
    { low<-mu[j]-4*s[j]
      upp<-mu[j]+4*s[j]
      y.out[,j]<-runif(nOut, min=low, max=upp)
    }
  }
  else { y.out<-NULL }
  return(list(y.out=y.out, nOut=nOut)) 
}

# Calculate the number of outliers
# numOutlier -- number/ratio of outliers. If numOutlier is a positive integer,
#   then numOutlier means the number of outliers. If numOutlier is a real
#   number between (0, 1), then numOutlier means the ratio of outliers, i.e. 
#   the number of outliers is equal to round(numOutlier*N), where N is the 
#   total number of non-outliers. If numOutlier is a real number greater than
#   1, then we round numOutlier to an integer. 
# N -- the total number of non-outliers
outlierSize<-function(numOutlier, N)
{
  if(numOutlier<0)
  {
    stop("The argument 'numOutlier' should be non-negative!\n")
  }
  N<-as.integer(N)
  if(N<0 || !is.integer(N))
  { stop("The number 'N' of non-outliers should be positive integer!\n") }

  if(numOutlier>=1) { nOut=floor(numOutlier) }
  else if(numOutlier>0) { nOut<-round(numOutlier*N) }
  else { nOut<-0 }
  return(nOut)
}


# Output the number of data points, the number of non-noisy variables, the
# number of noisy variables, the labels of noisy variables, the mean
# vectors, covariance matrices, theoretical and empirical separation indices
# and projection directions, asymptotic confidence lower bounds, rotation
# matrix, etc.
# b -- indicates the file is the b-th replicate
# fileName -- the first part of the log file names
# alpha -- the tuning parameter in the separation index
# sepVal -- the minimum separation index set as a priori
# G -- the number of cluster
# size -- cluster sizes
# N -- the total number of clusters
# p -- the total number of variables
# numNoisy -- the number of noisy variables
# noisySet -- the labels of noisy variables
# nOut -- the number of outliers
# thetaMat -- the mean vectors before rotation and randomizing orders
# s -- the array of covariance matrices before rotation and randomizing columns
# muMat -- the mean vectors after rotation and randomizing columns
# SigmaArray -- the array of covariance matrices after rotation and randomizing
#        columns
# Q -- the rotation matrix
# sepValMat -- the theoretical separation index matrix
# Jhat2 -- the empirical separation index matrix
# myprojDir -- the theoretical projection directions
# empProjDir -- the empirical projection directions
# egvaluesMat -- 'numClust' by 'p' eigenvalue matrix. Each row contains the 
#       'p' eigenvalues of a cluster, where 'numClust' is the number
#       of clusters, 'p' is the number of variables
# quiet -- a flag to switch on/off the outputs of intermediate results.
#       The default value is 'TRUE'.
# outputEmpirical -- indicates if empirical projection directions
#       and separation indices should be output to files.
#       These inforamtion usually are useful to check the cluster
#       structures. Hence, by default, 'outputEmpirical=TRUE'.
outputLog<-function(b, fileName, alpha, sepVal, G, size,
  N, p, numNoisy, noisySet, nOut, thetaMat, s, muMat, SigmaArray, Q, 
  sepValMat, Jhat2, myprojDir, empProjDir, egvaluesMat, quiet=TRUE,
  outputEmpirical=TRUE) 
{ # output log information, e.g. mean vectors and covariance matrices.
  if(!quiet)
  { cat(" *** Step 2.11:  Output results ***\n") }
  datafileName<-paste(fileName, "_", b, ".dat", sep="")
  logfileName<-paste(fileName, "_", b, ".log", sep="")
  tem<-paste("\n Log information for the data set ", datafileName, " >>", sep="")
  write.table(tem,append=FALSE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  # output parameters
  outputLogPar(logfileName,G,size,N,p,numNoisy,noisySet,nOut,alpha,
               sepVal, egvaluesMat)
  # output Mean vectors and covariance matrices before rotation 
  outputLogMeanCov(logfileName, G, thetaMat, s, flag=1)
  # output rotation matrix
  outputLogQ(logfileName, Q)
  # output Mean vectors and covariance matrices after rotation 
  outputLogMeanCov(logfileName, G, muMat, SigmaArray, flag=2)
  # output theoretical Separation index matrix and projection directions
  outputLogSepProj(logfileName, G, sepValMat, myprojDir)
  if(outputEmpirical)
  {
    # output empirical Separation index matrix and projection directions
    outputLogSepProjData(logfileName,G,alpha,Jhat2,empProjDir)
  }
  tem<-paste("\n *********** end file ", datafileName, " **********\n", sep="")
}

# output parameters
# logfileName -- file name of the log file
# G -- the number of cluster
# size -- cluster sizes
# N -- the total number of clusters
# p -- the total number of variables
# numNoisy -- the number of noisy variables
# noisySet -- the labels of noisy variables
# nOut -- the number of outliers
# alpha -- the tuning parameter in the separation index
# sepVal -- the minimum separation index which is set a priori
# egvaluesMat -- 'numClust' by 'p' eigenvalue matrix. Each row contains the 
#       'p' eigenvalues of a cluster, where 'numClust' is the number
#       of clusters, 'p' is the number of variables
outputLogPar<-function(logfileName,G,size,N,p,numNoisy,noisySet,nOut,alpha,
  sepVal, egvaluesMat)
{ 
  if(!is.numeric(alpha)) { tt1<-alpha } 
  else { tt1<-round(alpha, 3) }
  if(!is.numeric(sepVal)) { tt2<-sepVal } 
  else { tt2<-round(sepVal, 3) }
  tem<-paste("\n alpha = ", tt1, " sepVal = ", tt2, sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n Number of clusters >> ", G, sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n cluster sizes >>")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  write.table(t(size),append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n Number of data points >> ", N, sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n Number of dimensions >> ", p, sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n Number of noisy dimensions >> ", numNoisy, sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n noisy variables >>")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  write.table(t(noisySet),append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n Number of outliers >> ", nOut, sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n initial eigenvalues for each cluster >> \n (rows correspond to clusters)", sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  if(!is.numeric(egvaluesMat)) { tt<-egvaluesMat } 
  else { tt<-round(egvaluesMat, 3) }
  write.table(tt,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
}

# output Mean vectors and covariance matrices before/after rotation 
# logfileName -- name of the logfile
# G -- the number of clusters
# thetaMat -- the cluster center matrix 
# s -- the array of covariance matrices 
# flag -- flag =1 means before rotation and randomization
#         flag =2 means after rotation and randomization
outputLogMeanCov<-function(logfileName, G, thetaMat, s, flag=1)
{ p<-nrow(s[,,1])
  egvaluesMat<-matrix(0,nrow=G, ncol=p)
  rownames(egvaluesMat)<-paste("cluster", 1:G, sep="")
  colnames(egvaluesMat)<-paste("variable", 1:p, sep="")
  for(i in 1:G)
  { if(flag==1) # before rotation
    { tem<-paste("\n Mean vectors and covariance matrices before rotation --- ", i, "th cluster >>", sep="") }
    else { tem<-paste("\n Mean vectors and covariance matrices after rotation --- ", i, "th cluster >>", sep="") }
    write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
      row.names=FALSE,col.names=FALSE)
    tem<-paste("\n mean vector (mu'=Q*mu)>>")
    write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
      row.names=FALSE,col.names=FALSE)
    if(!is.numeric(thetaMat[i,]))
    { tt<-thetaMat[i,] 
      write.table(tt,append=TRUE, file=logfileName,quote=FALSE,
        row.names=FALSE,col.names=FALSE)
    } 
    else  
    { tt<-round(thetaMat[i,],3) 
      write.table(t(tt),append=TRUE, file=logfileName,quote=FALSE,
        row.names=FALSE,col.names=FALSE)
    }
    tem<-paste("\n covariance matrix (Sigma'=Q*Sigma*t(Q))>>\n")
    write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
      row.names=FALSE,col.names=FALSE)
    if(!is.numeric(s[,,i])) { tt<-s[,,i] } 
    else { tt<-round(s[,,i], 3) }
    write.table(tt,append=TRUE, file=logfileName,quote=FALSE,
      row.names=FALSE,col.names=FALSE)
    egvaluesMat[i,]<-eigen(s[,,i])$values
  }
  tem<-paste("\n Final eigenvalues for each cluster >> \n (rows correspond to clusters)", sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  if(!is.numeric(egvaluesMat)) { tt<-egvaluesMat } 
  else { tt<-round(egvaluesMat, 3) }
  write.table(tt,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
}

# output rotation matrix
# logfileName -- name of the logfile
# Q -- the rotation matrix
outputLogQ<-function(logfileName, Q)
{ tem<-paste("\n *************************************************\n", sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste(" ********************* rotation ******************\n", sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste(" *************************************************\n", sep="")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n The orthogonal matrix Q >>")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  if(!is.numeric(Q)) { tt<-Q } else { tt<-round(Q, 3) }
  write.table(tt,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
}

# output Separation index matrix and projection directions
# logfileName -- the name of the logfile
# G -- the number of clusters
# sepValMat -- the theoretical separation index matrix
# myprojDir -- the theoretical projection direction
outputLogSepProj<-function(logfileName, G, sepValMat, myprojDir)
{ tem<-paste("\n theoretical Separation index matrix >>\n")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  if(!is.numeric(sepValMat)) { tt<-sepValMat } 
  else { tt<-round(sepValMat, 3) }
  write.table(tt,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n theoretical projection directions >>")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  for(i in 1:(G-1))
  { for(j in (i+1):G)
    { tem<-paste("\n**** cluster ",i," and ", j,sep="");
      write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
        row.names=FALSE,col.names=FALSE)
      if(!is.numeric(as.vector(myprojDir[i,j,])))
      { tt<-as.vector(myprojDir[i,j,]) 
        write.table(tt,append=TRUE, file=logfileName,
          quote=FALSE, row.names=FALSE,col.names=FALSE)

      } 
      else  
      { tt<-round(as.vector(myprojDir[i,j,]),3) 
        write.table(t(tt),append=TRUE, file=logfileName,
          quote=FALSE, row.names=FALSE,col.names=FALSE)
      }
    }
  }
}

# output empirical Separation index matrix and projection directions
# logfileName -- the name of the logfile
# G -- the number of clusters
# alpha -- the tuning parameter in the separation index
# Jhat2 -- the empirical separation index matrix
# empProjDir -- the empirical projection directions
outputLogSepProjData<-function(logfileName,G,alpha, Jhat2,empProjDir)
{ tem<-paste("\n empirical Separation index matrix (alpha=",alpha,")>>\n")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  if(!is.numeric(Jhat2)) { tt<-Jhat2 } else { tt<-round(Jhat2, 3) }
  write.table(tt,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  tem<-paste("\n empirical projection directions >>")
  write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
    row.names=FALSE,col.names=FALSE)
  for(i in 1:(G-1))
  { for(j in (i+1):G)
    { tem<-paste("\n**** cluster ",i," and ", j,sep="");
      write.table(tem,append=TRUE, file=logfileName,quote=FALSE,
        row.names=FALSE,col.names=FALSE)
      if(!is.numeric(as.vector(empProjDir[i,j,])))
      { tt<-as.vector(empProjDir[i,j,])  
        write.table(tt,append=TRUE, file=logfileName,
          quote=FALSE, row.names=FALSE,col.names=FALSE)
      }
      else 
      { tt<-round(as.vector(empProjDir[i,j,]), 3)
        write.table(t(tt),append=TRUE, file=logfileName,
          quote=FALSE, row.names=FALSE,col.names=FALSE)
      }
    }
  }
}

# output data set
# b -- indicates the data set is the b-th replicate
# fileName -- the first part of the name of the data set
# y -- Nxp data matrix
# p -- the total number of variables
outputData<-function(b, fileName, y, p, outputDatFlag=TRUE)
{ tem<-paste("x", 1:p, sep="")
  if(outputDatFlag)
  { datafileName<-paste(fileName, "_", b, ".dat", sep="")
    write.table(y,file=datafileName,quote=FALSE,row.names=FALSE,col.names=tem)
  }
  colnames(y)<-tem
  return(y)
}


# Generate a Positive Definite Matrix 
# covMethod -- methods to generate random covariance matrices. 
#      'eigen' method first generates the eigenvalues of the covariance matrix,
#          then generate eigenvectors to construct the covariance matrix. 
#      'unifcorrmat' method first generates a correlation 
#          matrix via the method proposed by 
#           Joe H (2006). Generating random correlation matrices based on 
#           partial correlations. J. Mult. Anal. Vol. 97, 2177--2189
#          Then, it generate variances from the range 'rangeVar' to
#          construct the covariance matrix. 
#      There is no default method. 
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
# dim -- the dimension of this positive definite matrix
genPositiveDefMat<-function(dim, covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"), 
                            alphad=1, eta=1,
                            rangeVar=c(1,10), 
                            lambdaLow=1, ratioLambda=10)
{ 
  covMethod<-match.arg(arg=covMethod, choices=c("eigen", "onion", "c-vine", "unifcorrmat"))
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
  dim<-as.integer(dim)
  if(dim<1 || !is.integer(dim))
  { stop("The dimension 'dim' should be an integer greater than 1!\n") }

  if(covMethod=="eigen")
  { low<-lambdaLow
    upp<-lambdaLow*ratioLambda
  
    u<-matrix(0, dim,dim) # u is an diagonal matrix
    egvalues<-runif(dim,min=low,max=upp)
    diag(u)<-egvalues #the diagonal elements of u are positive
    Sigma<-u
    if(dim>1)
    { Q<-genOrthogonal(dim) # generate an orthogonal matrix 
      Sigma<-Q%*%u%*%t(Q) # the final positive definite matrix
    }
  } 
  else if (covMethod=="onion") {
    # use 'onion' method for random correlation matrix
    # and uniform distributions for random variances 
    # Reference: 
    # Joe H (2006). Generating random correlation matrices based on partial
    # correlations. J. Mult. Anal. Vol. 97, 2177--2189
    rr<-rcoronion(dim,eta)
    sigma2<-runif(dim, min=rangeVar[1], max=rangeVar[2])
    if(dim>1) { dd<-diag(sqrt(sigma2)) }
    else { dd<-sqrt(sigma2) }
    Sigma<-dd%*%rr%*%dd
    egvalues<-eigen(Sigma)$values

  }
  else if (covMethod=="c-vine"){
    # use 'c-vine' method for random correlation matrix
    # and uniform distributions for random variances 
    # Reference: 
    # Joe H (2006). Generating random correlation matrices based on partial
    # correlations. J. Mult. Anal. Vol. 97, 2177--2189
    rr<-rcorcvine(dim, eta)
    sigma2<-runif(dim, min=rangeVar[1], max=rangeVar[2])
    if(dim>1) { dd<-diag(sqrt(sigma2)) }
    else { dd<-sqrt(sigma2) }
    Sigma<-dd%*%rr%*%dd
    egvalues<-eigen(Sigma)$values
  }
  else 
  { # use hjrancor.R for random correlation matrix
    # and uniform distributions for random variances 
    # Reference: 
    # Joe H (2006). Generating random correlation matrices based on partial
    # correlations. J. Mult. Anal. Vol. 97, 2177--2189
    rr<-rcorrmatrix(dim, alphad=1)
    sigma2<-runif(dim, min=rangeVar[1], max=rangeVar[2])
    if(dim>1) { dd<-diag(sqrt(sigma2)) }
    else { dd<-sqrt(sigma2) }
    Sigma<-dd%*%rr%*%dd
    egvalues<-eigen(Sigma)$values
  }

  return(list(egvalues=egvalues, Sigma=Sigma))
}

# generate orthogonal matrix
# dim -- dimension
genOrthogonal<-function(dim)
{ 
  Q<-MOrthogonal(runif(dim))
  return(Q)
}

# Construct an orthogonal matrix whose first few columns are standardized 'M'
# where columns of 'M' are orthogonal.
# Here "standardized 'M'" means each its columns has length 1.
MOrthogonal<-function(M)
{
  # can set the parameter "tol" of "qr" to decide how small value should be 0
  tmp<-qr(M)
  Q<-qr.Q(tmp,complete=TRUE)
  if(is.vector(M)) { if(Q[1]*M[1]<0) Q<- -Q }
  else { if(Q[1,1]*M[1,1]<0) Q<- - Q }
  return(Q)
}



# Get a scalar 'A' such that the minimum separation index between clusters
# and their nearest neighboring clusters is equal to sepVal.
# The mean vector of the cluster 1 is (0,...,0). Other mean vectors
# have the form (0,...,0,A), (0,...,A,0),...,(A,0,...,0). 
# We want make sure the separation
# index between the cluster i and other clusters is at least sepVal,
# i=1,...,G, where G is the number of clusters.

# G -- the number of clusters
# sArray -- the array contains the covariance matrices of the clusters
# sepVal -- the minimum separation index which is set as a priori
# A2 -- the upper bound of A. 
# See documentation of genRandomClust for explanation of arguments:
# iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet
getA2<-function(G, sArray, sepVal, A2, 
                iniProjDirMethod=c("SL","naive"), 
                projDirMethod=c("newton", "fixedpoint"), 
                alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  G<-as.integer(G)
  if(G<1 || !is.integer(G))
  { stop("The number 'G' of clusters should be positive integer!\n")
  }
  if(sepVal<= -0.999 || sepVal>=0.999)
  { stop("The separation index 'sepVal' should be in (-0.999, 0.999)!\n")
  }
  if(A2<0)
  { stop("The initial upper bound 'A2' for the vertex-edge scalar 'A' should be positive!\n")
  }
  if(alpha<=0 || alpha>0.5)
  { stop("'alpha' should be between (0, 0.5]!\n")
  }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<1 || !is.integer(ITMAX))
  { stop("The maximum iteration allowed 'ITMAX' should be positive integer!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  # minimum scalar allowed so that the separation between a cluster and
  # its direct neighboring clusters is at least sepVal
  minA<- -3
  numNonNoisy<-ncol(sArray[,,1]) 
  # get the vertices of the simplex in numNonNoisy dimensional space
  # edge lengths of the simplex are all equal to 2
  vertex<-genShiftedVertexes(G, numNonNoisy) 
  d<-as.matrix(dist(vertex))
  myS<-1:G
  while(length(myS)>1)
  { i<-myS[1]
    myS<-myS[-1]
    myD<-which(abs(d[i, myS]-2)<eps) # direct neigbhors of cluster i
    myD<-myS[myD] 
    si<-sArray[,,i]
    vertexi<-vertex[i,]
    # make sure cluster j in myD and cluster i has separation index >= sepVal
    for(j in myD)
    { sj<-sArray[,,j]
      vertexj<-vertex[j,]; 
      # get interval [A1, A2] so that J^*_{ij}(A1)<0 and J^*_{ij}(A2)>0
      tmp<-getAInterval(vertexi,vertexj,si,sj,A2,sepVal,
                        iniProjDirMethod, projDirMethod, alpha,ITMAX,
                        eps, quiet)
      A1<-tmp$A1
      A2<-tmp$A2
      # find the value of A so that the separation index J^*_{ij}=sepVal.
      tmpA<-getA(A1, A2, vertexi, vertexj, si, sj, sepVal, 
                 iniProjDirMethod, projDirMethod, alpha, ITMAX,
                 eps, quiet)
      if(minA<tmpA) { minA<-tmpA }
    }
  }
  return(list(minA=minA, vertex=vertex))
}


# get interval [A1, A2] so that J^*_{ij}(A1)<0 and J^*_{ij}(A2)>0
# vertexi -- vertex 1
# vertexj -- vertex 2
# si -- covariance matrix for vertex 1
# sj -- covariance matrix for vertex 2
# A2 -- the upper bound of A. 
# sepVal -- the desired separation index between cluster i and cluster j 
# See documentation of genRandomClust for explanation of arguments:
# iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet
getAInterval<-function(vertexi, vertexj, si, sj, A2, sepVal,
                       iniProjDirMethod, projDirMethod, alpha=0.05, ITMAX=20,
                       eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))

  if(A2<0)
  { stop("The initial upper bound 'A2' for the vertex-edge scalar 'A' should be positive!\n")
  }
  if(sepVal<= -0.999 || sepVal>=0.999)
  { stop("The separation index 'sepVal' should be in (-1, 1)!\n") }
  if(alpha<=0 || alpha>0.5)
  { stop("'alpha' should be between (0, 0.5]!\n") }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<1 || !is.integer(ITMAX))
  { stop("The maximum iteration allowed 'ITMAX' should be positive integer!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  tmpsi<-min(diag(si), na.rm=TRUE) 
  tmpsj<-min(diag(sj), na.rm=TRUE)
  A1<-min(tmpsi, tmpsj, na.rm=TRUE)
  # get A1 so that the separation index between cluster i and j <0.
  while(1)
  { tmp<-funSepVal(A1, vertexi, vertexj, si, sj, sepVal, 
                   iniProjDirMethod, projDirMethod, alpha, 
                   ITMAX, eps, quiet)
    if(tmp>0) { A1<-A1/2 } else { break }
  }

  # get A2 so that the separation index between cluster i and j  >0. 
  tmpsi<-max(diag(si), na.rm=TRUE)
  tmpsj<-max(diag(sj), na.rm=TRUE)
  mystep<-max(tmpsi, tmpsj, na.rm=TRUE)
  while(1)
  { tmp2<-funSepVal(A2, vertexi, vertexj, si, sj, sepVal, 
                    iniProjDirMethod, projDirMethod, alpha, 
                    ITMAX, eps, quiet)
    if(tmp2<0) { A2<-A2+mystep } else { break }
  }
  return(list(A1=A1, A2=A2))
}

# estimate the value of A so that the separation index between cluster 1 and cluster 2 is sepVal
# A1 -- lower bound of A
# A2 -- upper bound of A
# i -- indicate we are dealing with the separation index between cluster 1 and cluster i
# vertex1 -- vertex 1
# vertex2 -- vertex i
# Sigma1 -- covariance matrix of the cluster 1
# Sigma2 -- covariance matrix of the cluster i
# sepVal -- the desired separation index between cluster 1 and cluster i 
# See documentation of genRandomClust for explanation of arguments:
# iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet
getA<-function(A1, A2, vertex1, vertex2, Sigma1, Sigma2, 
               sepVal, iniProjDirMethod=c("SL", "naive"),
               projDirMethod=c("newton", "fixedpoint"), alpha=0.05, ITMAX=20,
               eps=1.0e-10, quiet=TRUE)
{ 

  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(A1<0)
  { stop("The initial lower bound 'A1' for the vertex-edge scalar 'A' should be positive!\n")
  }
  if(A2<0)
  { stop("The initial upper bound 'A2' for the vertex-edge scalar 'A' should be positive!\n")
  }
  if(A1>A2)
  { stop("'A1' should be less than 'A2'!\n") }
  if(sepVal<= -0.999 || sepVal>=0.999)
  { stop("The separation index 'sepVal' should be between (-1, 1)!\n") }
  if(alpha<=0 || alpha>0.5)
  { stop("'alpha' should be between (0, 0.5]!\n") }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<1 || !is.integer(ITMAX))
  { stop("The maximum iteration allowed 'ITMAX' should be positive integer!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be between (0, 0.01]!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }
  newfit<-0
  class(newfit)<-"try-error"
  
  newfit<-try(
         tmp<-uniroot(f=funSepVal, lower=A1, upper=A2, 
               vertex1=vertex1, vertex2=vertex2, Sigma1=Sigma1, 
               Sigma2=Sigma2, sepVal=sepVal, 
               iniProjDirMethod=iniProjDirMethod,
               projDirMethod=projDirMethod, alpha=alpha, 
               ITMAX=ITMAX, eps=eps, quiet=quiet)
  )
  if(sum(class(newfit)=="try-error", na.rm=TRUE))
  { warning("Could not find suitable upper bound of 'myc'!\n 'myc' is set to be 'A2'!\n")
    A <-A2
  } else {
    A<-tmp$root
  }
  return(A)

#  A<-tmp$root
# return(A)
}

# given mui, Sigma1, Sigma2, first to calculate the separation index.
# then get the difference between this separation index and the given
# separation index.
# Sigma1 -- the covariance matrix of the first cluster. The mean vector of
#           the first cluster is 0.
# Sigma2 -- the covariance matrix of the second cluster. The mean vector of
#           the second cluser is (0,...,0, A, 0, ..., 0), where A is the
#           ith element.
# A -- see Sigma2
# i -- see Sigma2
# sepVal -- the given separation index
# See documentation of genRandomClust for explanation of arguments:
# iniProjDirMethod, projDirMethod, alpha, ITMAX, eps, quiet
funSepVal<-function(A, vertex1, vertex2, 
                    Sigma1, Sigma2, sepVal, 
                    iniProjDirMethod=c("SL", "naive"), 
                    projDirMethod=c("newton", "fixedpoint"), 
                    alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(sepVal<= -0.999 || sepVal>=0.999)
  { stop("The separation index 'sepVal' should be in (-0.999, 0.999)!\n") }
  if(alpha<=0 || alpha>0.5)
  { stop("'alpha' should be in (0, 0.5]!\n") }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<1 || !is.integer(ITMAX))
  { stop("The maximum iteration allowed 'ITMAX' should be positive integer!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }

  p<-ncol(Sigma1)
  mu1<-A*vertex1
  mu2<-A*vertex2
  # get the initial projection direction 
  iniProjDir<-getIniProjDirTheory(mu1, Sigma1, mu2, Sigma2, 
                            iniProjDirMethod, eps)
  # get the projection direction
  tmpa<-projDirTheory(iniProjDir, mu1, Sigma1, mu2, Sigma2, 
                         projDirMethod, alpha, ITMAX, eps, quiet)
  a<-tmpa$projDir
  # get separation index
  tmpsepVal<-tmpa$sepVal
  res<-tmpsepVal-sepVal
  return(res)
}


