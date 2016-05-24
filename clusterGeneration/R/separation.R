#
# Get initial projection direction
# if iniProjDirMethod="SL", iniProjDir=(Sigma1+Sigma2)^{-1}(mu_2-mu_1)
# else if iniProjDirMethod="naive", iniProjDir=(mu_2-mu_1)
# otherwise iniProjDir is randomly generated such that iniProjDir^%(mu2-mu1)>0
#
# mu1, Sigma1 -- mean vector and covariance matrix for cluster 1
# mu2, Sigma2 -- mean vector and covariance matrix for cluster 2
# iniProjDirMethod -- projDirMethod used to construct initial projection direction
#      projDirMethod="SL" => iniProjDir<-(Sigma1+Sigma2)^{-1}(mu2-mu1)
#      projDirMethod="naive" => iniProjDir<-(mu2-mu1)
# eps -- a small positive number used to check if a quantity is equal to zero.
#      'eps' is mainly used to check if an iteration algorithm converges.
getIniProjDirTheory<-function(mu1, Sigma1, mu2, Sigma2, 
                       iniProjDirMethod=c("SL", "naive"), 
                       eps=1.0e-10, quiet=TRUE)
{
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  if(iniProjDirMethod=="SL")
  { # Su and Liu (1993) projection direction, JASA 1993, vol88 1350-1355
    tmp<-Sigma1+Sigma2

    if(abs(det(tmp))<eps)
    { 
      t1<-sum(abs(as.vector(Sigma1)))
      t2<-sum(abs(as.vector(Sigma2)))
      if(abs((t1+t2))<eps) # Sigma_1=Sigma_2=0
      { if(!quiet)         
        { cat("Warning: Both covariance matrices are zero matrix!\n Naive direction is used instead!\n")
        }
        iniProjDir<-mu2-mu1
      } 
      else 
      { if(!quiet)
        { 
          cat("Warning: Generalized inverse is used!\n") 
        }
        iniProjDir<-ginv(Sigma1+Sigma2)%*%(mu2-mu1)
      }
    } 
    else { iniProjDir<-solve(Sigma1+Sigma2)%*%(mu2-mu1) }
  } 
  else 
  { # naive method
    iniProjDir<-as.vector(mu2-mu1)
  } 
  iniProjDir<-as.vector(iniProjDir)
  tmp<-as.vector(sqrt(crossprod(iniProjDir)))
  if(abs(tmp)>eps) 
  { 
    iniProjDir<-iniProjDir/tmp 
  } 
  else 
  { if(!quiet)
    {
      cat("Warning: iniProjDir=0!\n") 
    }
  }
  return(iniProjDir)
}


# Projection direction via iteration formula derived by 
# iteratively solving the equation d J(projDir) / d projDir = 0
# projDir -- projection direction at the t-th step
# output projection direction at the (t+1)-th step
# See documentation of getIniProjDirTheory for explanation of arguments:
# mu1, Sigma1, mu2, Sigma2, eps
getProjDirIter<-function(projDir, mu1, Sigma1, mu2, Sigma2, 
                         eps=1.0e-10, quiet=TRUE)
{
  tmp1<-as.numeric((abs(t(projDir)%*%Sigma1%*%projDir)))
  tmp2<-as.numeric((abs(t(projDir)%*%Sigma2%*%projDir)))
  b1ProjDir<-sqrt(tmp1)
  b2ProjDir<-sqrt(tmp2)
  diff<-mu2-mu1
  if(abs(tmp1)<eps && abs(tmp2)<eps)
  { if(!quiet) 
    { 
      cat("Warning: Both covariance matrices have rank zero!\n")
    }
    tt<-abs(diff)
    pos<-which(tt==max(tt, na.rm=TRUE))
    pos<-pos[1]
    projDir2<-rep(0,length(mu1))
    if(abs(tt[pos])>eps) # two clusters are totally separated
    { 
      projDir2[pos]<-1
      projDir2<-sign(diff[pos])*projDir2
    } 
    res<-list(stop=TRUE, projDir=projDir2)
    return(res)
  }
  if(abs(b1ProjDir)<eps) { tmp<-Sigma2/b2ProjDir } 
  else if (abs(b2ProjDir)<eps) { tmp<-Sigma1/b1ProjDir } 
  else { tmp<-Sigma1/b1ProjDir+Sigma2/b2ProjDir }
  inv<-ginv(tmp) 
  projDir2<-(b1ProjDir+b2ProjDir)*inv%*%diff
  res<-list(stop=FALSE, projDir=projDir2)
  return(res)
}

# Iteratively calculating optimal projection direction
# iniProjDir -- initial projection direction
# ITMAX  -- maximum iteration allowed
# quiet -- a flag to switch on/off the outputs of intermediate results.
#      The default value is 'TRUE'.
# See documentation of getIniProjDirTheory for explanation of arguments:
# mu1, Sigma1, mu2, Sigma2, eps
optimProjDirIter<-function(iniProjDir, mu1, Sigma1, mu2, Sigma2, 
                           ITMAX=20, eps=1.0e-10, quiet=TRUE)
{
  diff<-1.0e+300
  loop<-0
  denom<-sqrt(sum(iniProjDir^2))
  if(abs(denom)<eps)
  { if(!quiet) 
    { 
      cat("Warning: Initial projection direction 'iniProjDir'=0!\n")
    }
    projDir<-mu2-mu1
    denom<-sqrt(sum(projDir^2))
    if(abs(denom)<eps)
    { if(!quiet)
      { 
        cat("Warning: mu1=mu2\n") 
        cat("Optimal projection direction will be set to be zero\n")
      }
      return(rep(0, length(projDir)))
    }
    return((mu2-mu1)/denom) 
  }
  iniProjDir<-as.vector(iniProjDir/sqrt(sum(iniProjDir^2)))
  projDirOld<-iniProjDir
  while(1)
  {
    tmp<-getProjDirIter(projDirOld, mu1, Sigma1, mu2, Sigma2, eps, quiet) 
    projDir<-tmp$projDir
    stop<-tmp$stop
    if(stop==TRUE || sum(abs(projDir))<eps) { break }
    projDir<-as.vector(projDir/sqrt(sum(projDir^2)))
    diff<-max(abs(projDir-projDirOld), na.rm=TRUE)
    loop<-loop+1
    if(diff<eps) { break }
    if(loop>ITMAX)
    { if(!quiet)
      { 
        cat("Warning: Iterations did not converge!\n") 
      }
      break 
    }
    projDirOld<-projDir
  }  
  if(!quiet)
  { cat("number of iterations=", loop, " diff=", diff, "\n")
    cat("standardized iniProjDir>>\n"); print(iniProjDir); cat("\n");
    cat("standardized projDir>>\n"); print(projDir); cat("\n");
  }
  return(projDir)
}

# Separation index given a projection direction 'projDir'
# alpha -- tuning parameter for separation index to indicating the percentage 
#      of data points to downweight. We set 'alpha=0.05' like we set
#      the significance level in hypothesis testing as 0.05.
# See documentation of getIniProjDirTheory for explanation of arguments:
# mu1, Sigma1, mu2, Sigma2, eps
sepIndexTheory<-function(projDir, mu1, Sigma1, mu2, Sigma2, 
                         alpha=0.05, eps=1.0e-10, quiet=TRUE)
{
  if(as.vector(crossprod(projDir, mu2-mu1))<0)
  { if(!quiet) 
    {
      cat("Warning: 'projDir*(mu2-mu1)<0'! '-projDir' will be used!\n") 
    }
    projDir<- - projDir
  }
  # standardize projDir
  denom<-sqrt(sum(projDir^2))
  if(abs(denom)<eps)
  { if(!quiet)
    {
      cat("Warning: Projection direction 'projDir'=0!\n")
    }
    return(-1) 
  }
  projDir<-projDir/sqrt(sum(projDir^2))
  tmp1<-as.numeric((abs(t(projDir)%*%Sigma1%*%projDir)))
  tmp2<-as.numeric((abs(t(projDir)%*%Sigma2%*%projDir)))
  b1ProjDir<-sqrt(tmp1)
  b2ProjDir<-sqrt(tmp2)
  diff<-mu2-mu1
  if(abs(tmp1)<eps && abs(tmp2)<eps)
  { if(!quiet) 
    { 
      cat("Warning: Both covariance matrices have rank zero!\n")
    }
    tt<-sum(projDir*diff)
    if(abs(tt)<eps)
    { if(!quiet)
      {
        cat("Warning: Two cluster centers are the same!\n") 
      }
      return(0) # two clusters are totally overlaped
    } 
    else 
    { return(1) # two clusters are totally separated
    }
  }
  za<-qnorm(1-alpha/2)
  part1<-sum(projDir*diff)
  part2<-za*(b1ProjDir+b2ProjDir)
  d<-(part1-part2)/(part1+part2)
  return(d)
}


# Calculate separation index in one dimensional space, 
# given projected means and variances.
# make sure mu2>mu1 before using this function
# mu1, tau1 -- mle of the mean and standard devitation of cluster 1
# mu2, tau2 -- mle of mean and standard devitation of cluster 2
# alpha -- tuning parameter
sepIndex<-function(mu1, tau1, mu2, tau2, alpha=0.05, eps=1.0e-10)
{ Za<-qnorm(1-alpha/2) 
  L1<-mu1-Za*tau1; U1<-mu1+Za*tau1;
  L2<-mu2-Za*tau2; U2<-mu2+Za*tau2;
  denom<-U2-L1
  if(abs(denom)<eps)
  { sepVal<- -1 }
  else
  {
    sepVal<-(L2-U1)/(U2-L1); 
  }
  intercept<-(U1+L2)/2;
  return(list(sepVal=sepVal,intercept=intercept,L1=L1,U1=U1,L2=L2,U2=U2))
}

# Calculate the value of the separation index (data version)
# given a projection direction 'projDir'
# y1 -- data for cluster 1
# y2 -- data for cluster 2
# See documentation of sepIndexTheory for explanation of arguments:
# alpha, eps
sepIndexData<-function(projDir, y1, y2, alpha=0.05, eps=1.0e-10, quiet=TRUE)
{
  if(!(is.numeric(y1) || is.matrix(y1)))
  {
    stop("The argument 'y1' should be a numeric vector or matrix!\n")
  }
  if(!(is.numeric(y2) || is.matrix(y2)))
  {
    stop("The argument 'y2' should be a numeric vector or matrix!\n")
  }
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }

  if(is.vector(y1))
  { len<-length(y1) 
    mu1<-y1
    Sigma1<-matrix(0, nrow=len, ncol=len)
  } 
  else 
  { mu1<-apply(y1, 2, mean, na.rm=TRUE)
    Sigma1<-cov(y1)
  }
  if(is.vector(y2))
  { len<-length(y2) 
    mu2<-y2
    Sigma2<-matrix(0, nrow=len, ncol=len)
  } 
  else 
  { mu2<-apply(y2, 2, mean, na.rm=TRUE)
    Sigma2<-cov(y2)
  }

  res<-sepIndexTheory(projDir, mu1, Sigma1, mu2, Sigma2, alpha, eps, quiet)
  return(res)
}

# Use multiple initial projection directions to get multiple (local) optimal 
# projection directions and corresponding separation indices. 
# We also calculate separation indices for these initial projection directions
# directly.
# Then we choose the maximum separation index and its corresponding projection
# direction as output. This function is used to handle the case where
# one of or both of covariance matrices are singular
# We consider 2+2*p initial projection directions. The first two projection
# directions are 'SL' direction (Sigma1+Sigma2)^{-1}(Mu2-Mu1) and 
# 'naive' direction (Mu2-Mu1), respectively.
# The remaining projection directions are the eigenvectors of Sigma1 and Sigma2
# If iniProjDir^T(Mu2-Mu1)<0, iniProjDir = - iniProjDir
# Projection direction will be standardized so that it's length is equal to 1.
#
# ITMAX  -- maximum iteration allowed
# quiet -- a flag to switch on/off the outputs of intermediate results.
#      The default value is 'TRUE'.
# See documentation of getIniProjDirTheory and sepIndexTheory for explanation 
#   of arguments:
# mu1, Sigma1, mu2, Sigma2, alpha, eps
optimProjDirIterMulti<-function(mu1, Sigma1, mu2, Sigma2, alpha=0.05, 
                         ITMAX=20, eps=1.0e-10, quiet=TRUE)
{
  p<-length(mu1)
  num<-2+2*p
  projDirMat<-matrix(0, nrow=num, ncol=p)
  # get initial projection directions
  # Su and Liu (1993)'s direction
  projDirMat[1,]<-getIniProjDirTheory(mu1, Sigma1, mu2, Sigma2, "SL", 
                                      eps, quiet)
  # naive direction
  projDirMat[2,]<-getIniProjDirTheory(mu1, Sigma1, mu2, Sigma2, "naive", 
                                      eps, quiet)
  start<-2+1
  end<-2+p
  # eigenvectors of Sigma1
  projDirMat[start:end,]<-t(eigen(Sigma1)$vectors)
  start<-end+1
  end<-end+p
  # eigenvectors of Sigma2
  projDirMat[start:end,]<-t(eigen(Sigma2)$vectors)

  # for each initial direction, we get separation index
  pos<-0
  maxSep<- -2
  for(i in 1:num)
  { 
    tmpProjDir<-projDirMat[i,]
    if(sum(tmpProjDir*(mu2-mu1))<0)
    { tmpProjDir<- -tmpProjDir }
    sepVal1<-sepIndexTheory(tmpProjDir, mu1, Sigma1, mu2, Sigma2, 
                            alpha, eps, quiet)
    tmpProjDir2<-optimProjDirIter(tmpProjDir, mu1, Sigma1, mu2, 
                            Sigma2, ITMAX, eps, quiet) 
    sepVal2<-sepIndexTheory(tmpProjDir2, mu1, Sigma1, mu2, Sigma2, 
                            alpha, eps, quiet)
    
    if(sepVal1>maxSep)
    { maxSep<-sepVal1 # record maximum separation index
      projDirOpt<-tmpProjDir
    }
    if(sepVal2>maxSep)
    { maxSep<-sepVal2 # record maximum separation index
      projDirOpt<-tmpProjDir2
    }
  }
  return(list(projDir=projDirOpt, sepVal=maxSep))
}


#################################################################


# Finding optimal projection direction
# created by Weiliang Qiu, Jan. 23, 2005
# Obtain optimal projection direction for two sets of data points
# y1 -- n1 x p matrix for cluster 1
# y2 -- n2 x p matrix for cluster 2
# iniProjDirMethod -- method used to construct initial projection direction
#     iniProjDirMethod="SL" => iniProjDir<-(Sigma1+Sigma2)^{-1}(mu2-mu1)
#     iniProjDirMethod="naive" => iniProjDir<-(mu2-mu1)
# projDirMethod -- takes values "fixedpoint" or "newton"
#      "fixedpoint" means that we get optimal projection direction
#        by solving the equation d J(projDir) / d projDir = 0, where
#        J(projDir) is the separation index with projection direction 'projDir'
#      "newton" method means that we get optimal projection driection
#        by the method proposed in the appendix of Qiu and Joe (2006) 
#        "Generation of random clusters with specified degree of separation",
#        Journal of Classification Vol 23(2), 315--334 
# ITMAX -- maximum iteration allowed
# eps -- threshold for convergence criterion. If |new-old|<eps, then stop 
#       iterating
# alpha -- tuning parameter for separation index to indicating the percentage 
#      of data points to downweight. We set 'alpha=0.05' like we set
#      the significance level in hypothesis testing as 0.05.
# quiet -- indicating if intermediate results should be output
projDirData<-function(y1, y2, 
                iniProjDirMethod=c("SL", "naive"), 
                projDirMethod=c("newton", "fixedpoint"), 
                alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  if(!(is.numeric(y1) || is.matrix(y1)))
  {
    stop("The argument 'y1' should be a numeric vector or matrix!\n")
  }
  if(!(is.numeric(y2) || is.matrix(y2)))
  {
    stop("The argument 'y2' should be a numeric vector or matrix!\n")
  }
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<=0 || !is.integer(ITMAX))
  {
    stop("The maximum iteration number allowed 'ITMAX' should be a positive integer!\n")
  }
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  if(is.vector(y1))
  { len<-length(y1) 
    mu1<-y1
    Sigma1<-matrix(0, nrow=len, ncol=len)
  } 
  else 
  { mu1<-apply(y1, 2, mean, na.rm=TRUE)
    Sigma1<-cov(y1)
  }
  if(is.vector(y2))
  { len<-length(y2) 
    mu2<-y2
    Sigma2<-matrix(0, nrow=len, ncol=len)
  } 
  else 
  { mu2<-apply(y2, 2, mean, na.rm=TRUE)
    Sigma2<-cov(y2)
  }

  # get initial projection direction 'a' so that 'a^T*a=1'
  iniProjDir<-getIniProjDirTheory(mu1, Sigma1, mu2, Sigma2, 
                                  iniProjDirMethod, eps, quiet)
  # get optimal projection direction
  res<-projDirTheory(iniProjDir=iniProjDir,  mu1=mu1, Sigma1=Sigma1, 
                        mu2=mu2, Sigma2=Sigma2, 
                        projDirMethod=projDirMethod, alpha=alpha, 
                        ITMAX=ITMAX, eps=eps, quiet=quiet) 
  return(res)
}

# Obtain optimal projection direction for two sets of data points
# iniProjDir -- initial projection direction such that 't(iniProjDir)*(mu2-mu1)=1'
# mu1 -- mean vector of group 1
# Sigma1 -- covariance matrix of group 1
# mu2 -- mean vector of group 2
# Sigma2 -- covariance matrix of group 2
# projDirMethod -- takes values "fixedpoint" or "newton"
#      "fixedpoint" means that we get optimal projection direction
#        by solving the equation d J(projDir) / d projDir = 0, where
#        J(projDir) is the separation index with projection direction 'projDir'
#      "newton" method means that we get optimal projection driection
#        by the method proposed in the appendix of Qiu and Joe (2006) 
#        "Generation of random clusters with specified degree of separation",
#        Journal of Classification Vol 23(2), 315--334
# ITMAX -- maximum iteration allowed
# eps -- threshold for convergence criterion. If |new-old|<eps, then stop 
#       iterating
# alpha -- tuning parameter for separation index to indicating the percentage 
#      of data points to downweight. We set 'alpha=0.05' like we set
#      the significance level in hypothesis testing as 0.05.
# quiet -- indicating if intermediate results should be output
projDirTheory<-function(iniProjDir, mu1, Sigma1, mu2, Sigma2, 
                        projDirMethod=c("newton", "fixedpoint"), 
                        alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  # get eigenvalues and eigenvectors of the matrix 'Sigma1'
  det1<-det(Sigma1)
  # get eigenvalues and eigenvectors of the matrix 'Sigma2'
  det2<-det(Sigma2)
  if(!quiet)
  { cat("determinant of Sigma1=", det1, "\n")
    cat("determinant of Sigma2=", det2, "\n")
  }
  # one of covariance matrix is singular or use fixedpoint method
  # or iniProjDir^T*(mu2-mu1)=0.
  # For newton method, we require that iniProjDir^T*(mu2-mu1)=1.
  tt<-as.vector(crossprod(iniProjDir, mu2-mu1))
  if(abs(det1)<eps || abs(det2)<eps || abs(tt)<eps 
     || projDirMethod=="fixedpoint") 
  { if(!quiet) 
    { cat("Warning: Initial projection direction might be changed by the\n")
      cat("function 'projDirTheory()' because some cluster covariance \n")
      cat("matrices might be singular, or the argument 'projDirMethod' is\n")
      cat("specified as 'fixedpoint'!\n")
    }
    # search optimal projection direction using multiple starting projection
    # directions
    res<-optimProjDirIterMulti(mu1, Sigma1, mu2, Sigma2, 
                               alpha, ITMAX, eps, quiet) 
    return(res)
  }

  # normalize iniProjDir so that iniProjDir^T*(mu2-mu1)=1
  iniProjDir<-iniProjDir/tt
  # find matrix 'Q1' such that 't(Q1)*Sigma1*Q1=I_p'
  Q1<-Q1Fun(Sigma1)
  # find matrix 'Q2' and the scalar 'c1' such that 
  # 't(Q2)*t(Q1)*(mu2-mu1)=c1*e1',
  # where 'e1=c(1,0,...,0)'.
  tmp<-Q2c1Fun(Q1, mu1, mu2)
  Q2<-tmp$Q2
  c1<-tmp$c1
  # calculate the matrix 'V=t(Q2)*t(Q1)*Sigma2*Q1*Q2'
  V<-VFun(Q1, Q2, Sigma2)

  # inverse of 'Q1'
  iQ1<-solve(Q1)
  # inverse of 'Q2'
  iQ2<-solve(Q2)
  # projDir=iniProjDir = Q1*Q2*bini/c1
  bini<-c1*iQ2%*%iQ1%*%iniProjDir
  bini<-as.vector(bini)
  yini<-bini[-1]

  res<-newtonRaphson(yini, V, ITMAX, eps, quiet)
  y<-res$y
  code<-res$code
  b<-c(1, y)
  a<-Q1%*%Q2%*%b/c1
  a<-as.vector(a)

  # normalized a
  denom<-max(abs(a), na.rm=TRUE)
  if(abs(denom)<eps)
  { if(!quiet)
    {
      cat("Warning: projDir=0!\n") 
    }
    anorm<-a
  } 
  else 
  { #anorm<-a/max(abs(a), na.rm=TRUE)
    anorm<-a/as.vector(sqrt(crossprod(a)))
  }
  sepVal<-sepIndexTheory(anorm, mu1, Sigma1, mu2, Sigma2, alpha, eps, quiet)
  if(!quiet)
  { cat("normalized projection direction>>\n"); print(anorm);
  }

  return(list(projDir=anorm, sepVal=sepVal))
}

# iniProjDir proportional to (mu2-mu1)
# or iniProjDir proportional to LDA direction
# or randomly generate an initial projection direction
# y1 -- data for cluster 1
# y2 -- data for cluster 2
# iniProjDirMethod -- method used to construct initial projection direction
#           iniProjDirMethod="SL" => iniProjDir<-(Sigma1+Sigma2)^{-1}(mu2-mu1)
#           iniProjDirMethod="naive" => iniProjDir<-(mu2-mu1)
# eps -- threshold for convergence criterion. 
getIniProjDirData<-function(y1, y2, 
                            iniProjDirMethod=c("SL", "naive"), 
                            eps=1.0e-10, quiet=TRUE)
{
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  if(!(is.numeric(y1) || is.matrix(y1)))
  {
    stop("The argument 'y1' should be a numeric vector or matrix!\n")
  }
  if(!(is.numeric(y2) || is.matrix(y2)))
  {
    stop("The argument 'y2' should be a numeric vector or matrix!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }

  if(is.vector(y1))
  { len<-length(y1) 
    mu1<-y1
    Sigma1<-matrix(0, nrow=len, ncol=len)
  } 
  else 
  { mu1<-apply(y1, 2, mean, na.rm=TRUE)
    Sigma1<-cov(y1)
  }
  if(is.vector(y2))
  { len<-length(y2) 
    mu2<-y2
    Sigma2<-matrix(0, nrow=len, ncol=len)
  } 
  else 
  { mu2<-apply(y2, 2, mean, na.rm=TRUE)
    Sigma2<-cov(y2)
  }

  iniProjDir<-getIniProjDirTheory(mu1, Sigma1, mu2, Sigma2, 
                                  iniProjDirMethod, eps, quiet)
  return(iniProjDir) 
}


# find matrix 'Q1' such that 't(Q1)*Sigma1*Q1=I_p'
# require that Sigma1 is positive definite.
#
# What if Sigma1 is semi-positive definite, not positive definite?
# i.e., what if there are zero-value eigenvalues of Sigma1?
Q1Fun<-function(Sigma1)
{
  p<-nrow(Sigma1)
  # obtain eigenvalues and eigenvectors of 'Sigma1'
  eg<-eigen(Sigma1)
  eu<-eg$values # eigenvalues 'lambda_i, i=1,...,p'
  ieu2<-matrix(0,nrow=p, ncol=p)
  # diagonal elements are equal to '1/sqrt(lambda_i)'
  diag(ieu2)<-1.0/sqrt(eu)
  # columns of 'et' correspond to eigenvectors of 'Sigma1'
  et<-eg$vectors
  # 'Sigma1=et*diag(eu)*t(et)'
  # 'Q1=et*diag(eu)^{-1/2}*t(et)'
  Q1<-et%*%ieu2%*%t(et)
  return(Q1)
}

# find matrix 'Q2' and the scalar 'c1' such that 
# 't(Q2)*t(Q1)*(mu2-mu1)=c1*e1',
# where 'e1'=c(1,0,...,0)'
# 
# Q1 -- matrix such that 't(Q1)*Sigma1*Q1=I_p',
#       where 'Sigma1' is the covariance matrix for group 1
# mu1 -- mean vector for group 1
# mu2 -- mean vector for group 2
Q2c1Fun<-function(Q1, mu1, mu2)
{
  theta<-as.vector(t(Q1)%*%(mu2-mu1))
  # Q2=(q_{21}, q_{22}, ..., q_{2p})
  # t(q_{21})*theta=c1
  # t(q_{2i})*theta=0, i=2,\ldots, p
  Q2<-MOrthogonal(theta) # Q2 is an orthogonal matrix
  c1<-as.vector(Q2[,1]%*%theta)
  return(list(Q2=Q2, c1=c1))
}

# calculate the matrix 'V=t(Q2)*t(Q1)*Sigma2*Q1*Q2'
VFun<-function(Q1, Q2, Sigma2)
{
  V<-t(Q2)%*%t(Q1)%*%Sigma2%*%Q1%*%Q2
  return(V)
}


# quadratic form 't(y)*A*y'
quadraticFun<-function(y, A)
{
  y<-as.vector(y)
  A<-as.matrix(A)
  res<-t(y)%*%A%*%y
  res<-as.vector(res)
  return(res)
}

# the function g1() = t(y)*y+1
# y -- a (p-1) by 1 vector
g1Fun<-function(y)
{
  res<-crossprod(y)+1
  return(as.vector(res))
}

# the function g2() = (y+V22^{-1}v21)^T*V22*(y+V22^{-1}v21)+c2
# where c2=v11-v21^T*V22^{-1}*v21.
# By simplification, we can get
# g2() = y^T*V22*y+2*y^T*v21+v11
g2Fun<-function(y, V)
{
  tmp<-V2Fun(V)
  V22<-tmp$V22
  v21<-tmp$v21 
  v11<-V[1,1]
  # part1 = y^T*V22*y
  part1<-quadraticFun(y, V22)
  # part2 = 2*y^T*v21
  part2<-2*as.vector(crossprod(y, v21))
  res<-part1+part2+v11
  return(as.vector(res))
}

# the function g() = sqrt(g1(y))+sqrt(g2(y))
gFun<-function(y, V)
{
  g1<-g1Fun(y)
  g2<-g2Fun(y, V)
  g<-sqrt(g1)+sqrt(g2)
  return(g)
}

# the 1st derivative of g() = y/sqrt(g1(y))+(V22*y+v21)/sqrt(g2(y))
d1gFun<-function(y, V)
{
  g1<-g1Fun(y)
  g2<-g2Fun(y, V)
  tmp<-V2Fun(V)
  V22<-tmp$V22
  v21<-tmp$v21

  y<-as.vector(y)
  part1<-as.vector(y/sqrt(g1))
  part2<-(V22%*%y+v21)/sqrt(as.vector(g2))
  part2<-as.vector(part2)
  res<-part1+part2
  return(res)
}

# the 2nd derivative of g() = I/sqrt(g1(y))-y*y^T/[g1(y)*sqrt(g1(y))]
# + V22/sqrt(g2(y))-(V22*y+v21)*(V22*y+v21)^T/[g2(y)*sqrt(g2(y))]
d2gFun<-function(y, V)
{
  tmp<-V2Fun(V)
  V22<-tmp$V22
  v21<-tmp$v21
  g1<-as.vector(g1Fun(y))
  g2<-as.vector(g2Fun(y, V))
  p<-nrow(V)
  myI<-diag(p-1)
  part1<-(g1*myI-y%*%t(y))/as.vector(g1^(3/2))
  tmp<-V22%*%y+v21
  tmp<-tmp%*%t(tmp)
  part2<-(g2*V22-tmp)/as.vector(g2^(3/2))
  res<-part1+part2
  return(res)
}


# obtain V_{22} and v_{21}
V2Fun<-function(V)
{
  V22<-V[-1,-1]
  v21<-as.vector(V[,1])
  v21<-v21[-1]
  return(list(V22=V22, v21=v21))
}


# Solve projection direction via modified Newton-Raphson
# Convergence to eps=1.e-10 is  often in 3 iterations
# yini -- initial point for iterations
# V -- 'V=t(Q2)*t(Q1)*Sigma2*Q1*Q2'
# ITMAX -- the maximum iteration allowed
# eps -- convergence tolerance
# quiet -- a flag to switch on/off the outputs of intermediate results
newtonRaphson<-function(yini, V, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{
  code<-0
  loop<-0
  y<-yini
  while(1)
  {
    loop<-loop+1
    if(loop>ITMAX)
    { code<-1
      break
    }
    d1g<-d1gFun(y, V)
    d2g<-d2gFun(y, V)
    tem<-solve(d2g,d1g)
    tem<-as.vector(tem)
    ynew <- (y-tem)
    diff<-as.vector(sqrt(crossprod(tem)))
    if(diff<eps)
    { break }
    y<-ynew
    # modifed Newton-Raphson step, sometimes diverges without this
    while(diff>0.5) { tem<-tem/2; y<-y+tem; diff<-diff/2 }
  }
  if(code==1)
  { cat("warning! newtonRaphson() did not converge!\n") 
    cat("loop=", loop, " ITMAX=", ITMAX, " diff=", diff, 
      " eps=", eps, " code=", code, "\n")
  }

  if(!quiet)
  { cat("loop=", loop, " ITMAX=", ITMAX, " diff=", diff, 
      " eps=", eps, " code=", code, "\n")
  }
  return(list(y=y, code=code))
}

# Calculate the theoretical separation index and projection direction
# muMat -- cluster center matrix. muMat[i,] is the cluster center
#      for the i-th cluster
# SigmaArray -- array of covariance matrices. SigmaArray[,,i] is the 
#      covariance matrix for the i-th cluster
# iniProjDirMethod -- indicating the method to get initial projection direction 
#      By default, the sample version of the Su and Liu (SL) projection 
#      direction ((n1-1)*S1+(n2-1)*S2)^{-1}(mu2-mu1) is used,
#      where mui and Si are the mean vector and covariance
#      matrix of cluster i. Alternatively, the naive projection
#      direction (mu2-mu1) is used. 
#      Su and Liu (1993) JASA 1993, vol. 88 1350-1355.
# projDirMethod -- takes values "fixedpoint" or "newton"
#      "fixedpoint" means that we get optimal projection direction
#        by solving the equation d J(progDir) / d progDir = 0, where
#        J(progDir) is the separation index with projection direction 'progDir'
#      "newton" method means that we get optimal projection driection
#        by the method proposed in the appendix of Qiu and Joe (2006) 
#        "Generation of random clusters with specified degree of separation",
#        Journal of Classification Vol 23(2) 315--334
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
getSepProjTheory<-function(muMat, SigmaArray, 
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
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<=0 || !is.integer(ITMAX))
  {
    stop("The maximum iteration number allowed 'ITMAX' should be a positive integer!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  if(!quiet)
  { cat(" *** Step 2.10.1:  Calculate the theoretical separation index matrix and projection directions.***\n") }
  G<-nrow(muMat) # number of clusters
  p<-ncol(muMat) # number of dimensions
  projDirArray<-array(0, c(G,G,p))
  dimnames(projDirArray)<-list(NULL, NULL, paste("variable", 1:p, sep=""))
  sepValMat<-matrix(-1, nrow=G, ncol=G)
  rownames(sepValMat)<-paste("cluster", 1:G, sep="")
  colnames(sepValMat)<-paste("cluster", 1:G, sep="")
  for(i in 2:G)
  { mui<-muMat[i,] 
    si<-SigmaArray[,,i]
    for(j in 1:(i-1))
    { muj<-muMat[j,]
      sj<-SigmaArray[,,j]

      #**** begin 1 WQ 09/22/2007
      # in case some variances are equal to zero
      dsi<-as.numeric(diag(si))
      dsj<-as.numeric(diag(sj))
      tmppos<-which(abs(dsi-dsj)<eps)
      if(length(tmppos)>0)
      {
        # if difference of means is also equal to zero, then this dimension is noisy
        # we should set the corresonding elements in the projection direction as zero
        myset<-1:p
        # variables have zero variances in both clusters
        myset2<-myset[tmppos]
        diffmu<-abs(mui-muj)

        diffmu2<-diffmu[myset2]
        tmppos2<-which(diffmu2<eps)

        tmpn<-length(diffmu2)
        tmplen<-length(tmppos2)
        if(tmplen>0 && tmplen==tmpn) 
        {
          if(tmpn==p)
          { # all varialbes are noisy, i.e., the points in the two clusters have
            # the same coordinates.
            # then any projection direction will produce the same results.
            # So we set the optimal direction as (1,0,0,0...0)
            a<-rep(0, p)
            a[1]<-1
            tmpaSepVal<- -1
          } else if(tmpn==(p-1)) {
            # only one variable is non-noisy
            myset3<-myset[-tmppos]
            a<-rep(0, p)
            a[myset3]<-1
            nui<-mui[myset3]
            nuj<-muj[myset3]
            taui<-si[myset3, myset3]
            tauj<-sj[myset3, myset3]
            tmpaSepVal<-sepIndex(nui, taui, nuj, tauj, alpha, eps)
          } else { # we get projection direction for other dimensions
            myset3<-myset[-tmppos] 
            nui<-mui[myset3]
            nuj<-muj[myset3]
            taui<-si[myset3, myset3]
            tauj<-sj[myset3, myset3]
          
            # get the initial projection direction
            iniProjDir<-getIniProjDirTheory(nui, taui, nuj, tauj, 
                                            iniProjDirMethod, eps, quiet)
            # get the projection direction
            tmpa<-projDirTheory(iniProjDir,  nui, taui, nuj, tauj, 
                                projDirMethod, alpha, ITMAX, eps, quiet)
            a2<-tmpa$projDir
            tmpaSepVal<-tmpa$sepVal
          
            a<-rep(0, p)
            a[myset3]<-a2
          }
        }  else { # some variables are non-noisy, some variable are noisy
          # we choose a dimension has the maximum separation
          tmppos3<-which(diffmu2==max(diffmu2, na.rm=TRUE))[1]
          mydim<-myset2[tmppos3]
          a<-rep(0, tmpn) 
          a[mydim]<-1
          tmpaSepVal<-1
        }
      } else { # In all variables, variance of two clusters are not all equal to zero

        # get the initial projection direction
        iniProjDir<-getIniProjDirTheory(mui, si, muj, sj, 
                                        iniProjDirMethod, eps, quiet)
        # get the projection direction
        tmpa<-projDirTheory(iniProjDir,  mui, si, muj, sj, 
                            projDirMethod, alpha, ITMAX, eps, quiet)
        a<-tmpa$projDir
        tmpaSepVal<-tmpa$sepVal
      }
      projDirArray[i,j,]<-a
      projDirArray[j,i,]<- -a
      # get the separation index
      sepValMat[i,j]<-tmpaSepVal
      sepValMat[j,i]<-tmpaSepVal
      #**** end 1 WQ 09/22/2007
    }
  }
  return(list(sepValMat=sepValMat, projDirArray=projDirArray))
}

# get empirical separation indices and projection directions
# y -- the Nxp data matrix
# cl -- the memberships of data points
# iniProjDirMethod -- indicating the method to get initial projection direction 
#      By default, the sample version of the Su and Liu (SL) projection 
#      direction ((n1-1)*S1+(n2-1)*S2)^{-1}(mu2-mu1) is used,
#      where mui and Si are the mean vector and covariance
#      matrix of cluster i. Alternatively, the naive projection
#      direction (mu2-mu1) is used. 
#      Su and Liu (1993) JASA 1993, vol. 88 1350-1355.
# projDirMethod -- takes values "fixedpoint" or "newton"
#      "fixedpoint" means that we get optimal projection direction
#        by solving the equation d J(progDir) / d progDir = 0, where
#        J(progDir) is the separation index with projection direction 'progDir'
#      "newton" method means that we get optimal projection driection
#        by the method proposed in the appendix of Qiu and Joe (2006) 
#        "Generation of random clusters with specified degree of separation",
#        Journal of Classification Vol 23(2) 315--334
# alpha -- tuning parameter for separation index to indicating the percentage 
#       of data points to downweight. We set 'alpha=0.05' like we set
#       the significance level in hypothesis testing as 0.05.
# ITMAX -- when calculating the projection direction, we need iterations. 
#      ITMAX gives the maximum iteration allowed.
#      The actual #iterations is usually much less than the default value 20.
# eps -- convergence tolerance
# quiet -- a flag to switch on/off the outputs of intermediate results.
#      The default value is 'TRUE'.
getSepProjData<-function(y, cl, 
                         iniProjDirMethod=c("SL", "naive"), 
                         projDirMethod=c("newton", "fixedpoint"), 
                         alpha=0.05, ITMAX=20, eps=1.0e-10, quiet=TRUE)
{ 
  # QWL: should handle univariate case in the next version!
  if(!(is.matrix(y) || is.data.frame(y)))
  { stop("Error! The argument 'y' should be a matrix or data frame with row number (observation number) greater than 1!\n")
  }
  if(nrow(y)!=length(cl))
  { stop("The row number of 'y' does not match the length of 'cl'!\n") }
  p<-ncol(y)
  if(alpha<=0 || alpha>0.5)
  {
    stop("The tuning parameter 'alpha' should be in the range (0, 0.5]!\n")
  }
  iniProjDirMethod<-match.arg(arg=iniProjDirMethod, choices=c("SL", "naive"))
  projDirMethod<-match.arg(arg=projDirMethod, choices=c("newton", "fixedpoint"))
  ITMAX<-as.integer(ITMAX)
  if(ITMAX<=0 || !is.integer(ITMAX))
  {
    stop("The maximum iteration number allowed 'ITMAX' should be a positive integer!\n")
  }
  if(eps<=0 || eps > 0.01)
  {
    stop("The convergence threshold 'eps' should be in (0, 0.01]!\n")
  }
  if(!is.logical(quiet))
  {
    stop("The value of the quiet indicator 'quiet' should be logical, i.e., either 'TRUE' or 'FALSE'!\n")
  }

  if(!quiet)
  { cat(" *** Step 2.10.2:  Calculate the empirical separation index matrix and projection directions.***\n") }
  u.cl<-sort(unique(cl))
  G<-length(u.cl)
  if(G<2)
  { stop("Error! There is only one cluster in the data set!\n") }
  # calculate the matrix of mean vectors and the array of covariance matrices
  muMat<-matrix(0, nrow=G, ncol=p)
  SigmaArray<-array(0, c(p,p,G)) 
  dimnames(SigmaArray)<-list(NULL, NULL, paste("cluster", 1:G, sep=""))
  for(i in 1:G)
  { yi<-y[which(cl==u.cl[i]),,drop=FALSE]
    muMat[i,]<-apply(yi, 2, mean, na.rm=TRUE)
    if(nrow(yi)>1) { SigmaArray[,,i]<-cov(yi) } 
    else { SigmaArray[,,i]<-matrix(0, nrow=p, ncol=p) }
  }
  res<-getSepProjTheory(muMat, SigmaArray, 
                        iniProjDirMethod, projDirMethod,
                        alpha, ITMAX, eps, quiet)
  return(res)
}

