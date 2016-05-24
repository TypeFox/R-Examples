#===========================================================================#
# caTools - R library                                                       #
# Copyright (C) 2005 Jarek Tuszynski                                        #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

LogitBoost = function(xlearn, ylearn, nIter=ncol(xlearn))
#   An implementation of the LogitBoost classification algorithm with
#   decision stumps as weak learners. 
#   
# Input:
#   xlearn - A dataset, whose n rows contain the training instances.
#   ylearn - Class labels of dataset
#   nIter - An integer, describing the number of iterations for
#     which boosting should be run. 
#     
# Output:
#   Stump - list of decision stumps used:
#            column 1: contains feature numbers or each stump
#            column 2: contains threshold
#            column 3: contains bigger/smaller info: 1 means (col>thresh ? class_1 : class_2);  -1 means (col>thresh ? class_2 : class_1)
#
# Writen by Jarek Tuszynski - SAIC (jaroslaw.w.tuszynski@saic.com)
# This code was addapted from logitboost.R function written by Marcel Dettling
# See "Boosting for Tumor Classification of Gene Expression Data", Dettling and Buhlmann (2002), 
#   available on the web page http://stat.ethz.ch/~dettling/boosting.html
{
  lablist = sort(unique(ylearn))     # different classes in label array
  nClass  = length(lablist)          # number of different classes in label array
  
  if (nClass>2) {                    # Multi class version uses 2-class code ...
    Stump = NULL                     # ... recursivly
    for (jClass in 1:nClass) {
      y = as.numeric(ylearn!=lablist[jClass]) # lablist[jClass]->1; rest->0
      Stump = cbind(Stump, LogitBoost(xlearn, y, nIter)$Stump)
    }
    object = list(Stump=Stump, lablist=lablist) # create LogitBoost object
    class(object) <- "LogitBoost"
    return(object)
  }
  
  Mask = is.na(xlearn)                # any NA in test data will ... 
  if (any(Mask)) xlearn[Mask] = Inf   # ... be changed to +Inf
  ylearn = as.numeric(ylearn!=lablist[1]) # change class labels to boolean 
  nLearn = nrow(xlearn)               # Length of training data         
  nFeat  = ncol(xlearn)               # number of features to choose from         
 
  # Array Initialization
  f      = 0                          # range -1 or 1          
  p      = numeric(nLearn)+1/2        # range [0,1]
  Stump  = matrix(0, nIter,3)         # will hold the results
  colnames(Stump) = c("feature", "threshhold", "sign")
  Thresh =            matrix(0, nLearn, nFeat)  # sorted xlearn each fearure at a time
  Index  = matrix(as.integer(0), nLearn, nFeat) # order of samples in Thresh
  Mask   = matrix(as.logical(0), nLearn, nFeat) # mask of unique samples in Thresh 
  repts  = as.logical(numeric(nFeat))           # for each feature: are all samples unique
 
  # sort all columns and store them in Thresh; Index stores original positions
  for (iFeat in 1:nFeat) {  
    x = sort(xlearn[,iFeat], index=TRUE)
    Thresh[,iFeat] = x[[1]]
    Index [,iFeat] = x[[2]]
    Mask  [,iFeat] = c((diff(Thresh[,iFeat])!=0), TRUE) 
    repts [ iFeat] = !all(Mask[,iFeat])
  }

  # Boosting Iterations
  jFeat = 0
  for (m in 1:nIter) {
    # Computation of working response and weights
    w = pmax(p*(1-p), 1e-24) # weight for each sample
    z = (ylearn-p)/w
    w = w/sum(w)             # normalize to 1

    # older version similar to rpart (more readable but slower) left as documentation
      # MinLS = 1000                  # initialize search for minimum Least-Square
      # for (iFeat in 1:nFeat) {      # for each feature do: for each spliting point calculate least square regresion of z to xlearn with weights w.
        # Col = xlearn[,iFeat]        # sort all columns and store each one in thr
        # thr = Thresh[,iFeat]        # sort all columns and store each one in thr
        # for (j in 1:nLearn) {       # for every possible threshold
          # ff = 2*(Col<=thr[j])-1     # for every point in the column if (col<=Thresh) then f=1 else f=-1 
          # LS1[j] = sum(w*(z-ff)^2)   # Least Square array for every possible theshold  ( f = (col<=Thresh ? 1 : -1) )
          # LS2[j] = sum(w*(z+ff)^2)   # Least Square array for every possible theshold  after swaping output classes ( f = (col<sp ? -1 : 1) )
        # }
        # iLS1 = which.min(LS1)
        # iLS2 = which.min(LS2)
        # vLS1 = LS1[iLS1]
        # vLS2 = LS2[iLS2]      
        # if (MinLS>vLS1) { stump=c(iFeat, thr[iLS1],  1); MinLS=vLS1; }
        # if (MinLS>vLS2) { stump=c(iFeat, thr[iLS2], -1); MinLS=vLS2; }
      # }
      
    # Same as code above but faster:
    # for each feature do: for each spliting point calculate least square regresion of z to xlearn with weights w.
    # This part of the code is my version of rpart object.
    # LS1 is least square array LS1 = sum((w.*(z-f)).^2) where f = (col<=sp ? 1 : -1) LS1 is calculated for all spliting points
    # LS2 is least square array LS2 = sum((w.*(z-f)).^2) where f = (col<=sp ? -1 : 1) LS2 is calculated for all spliting points
    # for each spliting point col(idx) we will take one sample and change its sign, that will change current least square:
    # LS(i) = LS(i) - ( w(idx(i)) .* ( z(idx(i)) - (-1) ) ).^2   +   ( w(idx(i)) .* ( z(idx(i)) - 1 ) ).^2 what can be simplified to:
    # LS(i) = LS(i) - 4*(w(idx(i)).^2) .*  z(idx(i))
    ls1 = sum(w*(z+1)^2)     # Least Square value for left-most theshold  ( f = -1 )
    ls2 = sum(w*(z-1)^2)     # Least Square value for left-most theshold after swaping output classes ( f =  1 )
    MinLS = max(ls1,ls2)     # initialize search for minimum Least-Square
    wz  = 4 * w * z          # precompute vector to be used later
    for (iFeat in 1:nFeat) { # for each feature do: for each spliting point calculate least square regresion of z to xlearn with weights w.
      if (iFeat==jFeat) next # Prevent the simplest cycle
      Col = Thresh[,iFeat]   # get one column of sorted data
      LS  = cumsum(wz[Index[,iFeat]])  # find offset to Least Square value for every possible theshold
      if (repts[iFeat]) {    # if any repeating thresholds than delete all but last LS value
        mask = Mask[,iFeat]  # use mask of unique samples 
        Col = Col[mask]      # delete non-unique thresholds since they can cause errors
        LS  = LS [mask]      # delete coressponding Least Square values         
      }
      iLS1 = which.max(LS)   # min of LS1=ls1-LS - Least Square value for every possible theshold  ( f = (col<=Thresh ? 1 : -1) )
      iLS2 = which.min(LS)   # min of LS2=ls2+LS - Least Square value for every possible theshold after swaping output classes ( f = (col<Thresh ? -1 : 1) )
      vLS1 = ls1-LS[iLS1]
      vLS2 = ls2+LS[iLS2]
      if (MinLS>vLS1) { stump=c(iFeat, Col[iLS1],  1); MinLS=vLS1; }
      if (MinLS>vLS2) { stump=c(iFeat, Col[iLS2], -1); MinLS=vLS2; }
    }
 
    # =================
    # Fitting the tree
    # =================
    Stump[m,] = stump                                     # if stump[3]>0  f(i) += (xlearn(i,stump[1])<=stump[2] ? 1 : -1)
    jFeat = stump[1]
    f = f + stump[3] * (2*( xlearn[,jFeat]<=stump[2] )-1) # if stump[3]<0  f(i) += (xlearn(i,stump[1])> stump[2] ? 1 : -1)
    p = 1/(1+exp(-f))    # Updating and probabilities - range (0, 1]
    
    y = (f>0)
    y[f==0] = 0.5
    Conv = sum(abs(ylearn-y)) # keep track of error rate
  }
  object = list(Stump=Stump, lablist=lablist)  # create LogitBoost object
  class(object) <- "LogitBoost"
  return(object)
}


predict.LogitBoost = function(object, xtest, type = c("class", "raw"), nIter=NA, ...)
#   An implementation of the LogitBoost classification algorithm with
#   decision stumps as weak learners. 
#   
# Input:
#   xtest - A dataset, whose n rows contain samples to be classified.
#   ytest - Class labels of dataset (if given than Conv will return error rates as function of iterations)
#     
# input:
#   Stump - list of decision stumps used:
#            column 1: contains feature numbers or each stump
#            column 2: contains threshold
#            column 3: contains bigger/smaller info: 1 means (col>Thresh ? class_1 : class_2);  -1 means (col>Thresh ? class_2 : class_1)
#
# Writen by Jarek Tuszynski - SAIC (jaroslaw.w.tuszynski@saic.com)
# This code was addapted from logitboost.R function written by Marcel Dettling
# See "Boosting for Tumor Classification of Gene Expression Data", Dettling and Buhlmann (2002), 
#   available on the web page http://stat.ethz.ch/~dettling/boosting.html
{
  type    = match.arg(type)
  Stump   = object$Stump
  lablist = object$lablist
  if (is.na(nIter)) nIter = nrow(Stump)
  else nIter  =  min(nIter, nrow(Stump))
  nTest   = nrow(xtest)
  nClass  = ncol(Stump)/3
  if (nClass==1) nClass=2
  Prob    = matrix(0,nTest, nClass)
  colnames(Prob) = lablist
  if (nClass>2) {                     # multi class problem
    object$lablist = c(1,2)           # generic labels 
    for(iClass in 1:nClass) {
      object$Stump = Stump[,3*iClass + (-2:0)]
      prob = predict.LogitBoost(object, xtest, type="raw", nIter=nIter)
      Prob[,iClass] = prob[,1] # probability that "sample belongs to iClass" is true 
    }
  } else {                            # two class problem
    Mask    = is.na(xtest)            # any NA in test data will... be changed 
    if (any(Mask)) xtest[Mask] = Inf  # be changed to +Inf
    f = numeric(nTest)           
    for (iter in 1:nIter) {
      iFeat  = Stump[iter,1]
      thresh = Stump[iter,2]
      Sign   = Stump[iter,3]
      f = f + Sign*(2*(xtest[,iFeat]<=thresh)-1)
    }
    Ytest = (f>0)
    Ytest[f==0] = 0.5
    prob = 1/(1+exp(-f))
    Prob[,1] = 1-prob # probability that "sample belongs to 1" 
    Prob[,2] =   prob # probability that "sample belongs to 2" 
  }
  if (type=="raw") RET=Prob        # request to return raw probabilities
  else {                           # otherwise assign labels
    ord = t(apply(-Prob, 1, order))# find order of sorted Probs 
    RET = lablist[ord[,1]]         # find label with highest Prob
    Prob = -t(apply(-Prob,1,sort)) # sort Probs
    RET[Prob[,1]==Prob[,2]] = NA   # in case of ties return NA's
  }
  return (RET)
}

