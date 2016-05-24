#-------------------------
# Benjamin Risk
#Match based on L2 distances
matchICA=function(S,template,M=NULL) {
  n.comp=ncol(S)
  n.obs=nrow(S)
  if(n.comp>n.obs) warning('Input should be n x d')
  if(!all(dim(template)==dim(S))) warning('Template should be n x d')
  S = t(S)
  template = t(template)
  l2.mat1=matrix(NA,nrow=n.comp,ncol=n.comp)  
  l2.mat2=l2.mat1
  for (j in 1:n.comp) {
    for (i in 1:n.comp) {
      l2.mat1[i,j]=sum((template[i,]-S[j,])^2)/n.obs
      l2.mat2[i,j]=sum((template[i,]+S[j,])^2)/n.obs
    }
  }
  l2.mat1=sqrt(l2.mat1)
  l2.mat2=sqrt(l2.mat2)
  l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
  map=as.vector(solve_LSAP(l2.mat))
  l2.1=diag(l2.mat1[,map])
  l2.2=diag(l2.mat2[,map])
  sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2) 
  perm=diag(n.comp)[,map]%*%diag(sign.change)
  
  s.perm=t(perm)%*%S
  if(!is.null(M)) {
    M.perm=t(M)%*%perm
    return(list(S=t(s.perm),M=t(M.perm)))
  }  else {
    t(s.perm)
  }
}

#-------------------------------------
# Match mixing matrices:
# This function does not require M to be square:
# 26 September 2014 BRisk: Modified function for M1 d x p and M2 q x p where d may not be equal to q
frobICA<-function(M1=NULL,M2=NULL,S1=NULL,S2=NULL,standardize=FALSE) {
  #MODEL: X = S M + E, so M is d x p
  #standardize: if standardize==TRUE, then standardizes rows of M1 and M2 
  #to have unit norm; if using S1 and S2, standardizes columns to have unit variance.
  #standardize=TRUE makes the measure scale invariant.
  
  tfun = function(x) all(x==0)
  if(is.null(M1) && is.null(M2) && is.null(S1) && is.null(S2)) stop("need to supply either M1 and M2 or S1 and S2")
  if(!is.null(M1) && !is.null(M2) && !is.null(S1) && !is.null(S2)) {
    stop("provide either (M1 and M2) or (S1 and S2) but not both (M1,M2) and (S1,S2)")  
  }
  if(!is.null(M1) && nrow(M1) > ncol(M1)) stop("The input appears to be S1 and S2, but the arguments were not specified; re-run with S1=<object> and S2=<object>")
  
  if(is.null(M1)) {
    nS = nrow(S1)
    if(nS!=nrow(S2)) stop('S1 and S2 must have the same number of rows')
    if(sum(apply(S1,2,tfun)) + sum(apply(S2,2,tfun))) stop('frobICA not defined when S1 or S2 has a column of all zeros')
    if(standardize) {
      S1 = scale(S1)
      S2 = scale(S2)
    }
    p = ncol(S1)
    q = ncol(S2)
    if(p < q) {
      S1 = cbind(S1,matrix(0,nS,(q-p)))
    }
    if(q < p) {
      S2 = cbind(S2,matrix(0,nS,(p-q)))
    }
    Stemp = matchICA(S=S1,template=S2)
    n.comp = max(q,p)
    indices = c(1:n.comp)[!(apply(Stemp,2,tfun) | apply(S2,2,tfun))]
    return(sqrt(sum((Stemp[,indices] - S2[,indices])^2))/sqrt(nS*min(p,q)))
  } 
  
  else {
    if(sum(apply(M1,1,tfun)) + sum(apply(M2,1,tfun))) stop('frobICA not defined when M1 or M2 has a row of all zeros')
    if(standardize) {
      temp = diag((diag(M1%*%t(M1)))^(-1/2))
      M1 = temp%*%M1
      temp = diag((diag(M2%*%t(M2)))^(-1/2))
      M2 = temp%*%M2
    }
  p = ncol(M1)
  if(p!=ncol(M2)) stop("M1 and M2 must have the same number of columns")
  d = nrow(M1)
  q = nrow(M2)
  n.comp=max(d,q)
  if(n.comp > p) warning("M should be d x p")
  if(d<q) {
    M1 = rbind(M1,matrix(0,(q-d),p))
  }
  if(q<d) {
    M2 = rbind(M2,matrix(0,(d-q),p))
  }
  l2.mat1=l2.mat2=matrix(NA,nrow=n.comp,ncol=n.comp)  
  for (j in 1:n.comp) {
    for (i in 1:n.comp) {
      #since signs are arbitrary, take min of plus and minus:
      l2.mat1[i,j]=sum((M2[i,]-M1[j,])^2)
      l2.mat2[i,j]=sum((M2[i,]+M1[j,])^2)
    }
  }
  l2.mat1=sqrt(l2.mat1)
  l2.mat2=sqrt(l2.mat2)
  #take the min of plus/min l2 distances. This is okay because solve_LSAP is one to one
  l2.mat=l2.mat1*(l2.mat1<=l2.mat2)+l2.mat2*(l2.mat2<l2.mat1)
  map=as.vector(solve_LSAP(l2.mat))
  #retain relevant l2 distances:
  l2.1=diag(l2.mat1[,map])
  l2.2=diag(l2.mat2[,map])
  #sign.change is for re-ordered matrix 2
  sign.change=-1*(l2.2<l2.1)+1*(l2.1<=l2.2) 
  perm=diag(n.comp)[,map]%*%diag(sign.change)
  M.perm=t(perm)%*%M1
  indices = c(1:n.comp)[!(apply(M.perm,1,tfun) | apply(M2,1,tfun))]
  return(sqrt(sum((M.perm[indices,]-M2[indices,])^2))/sqrt(p*min(d,q)))
  } 
}


#----------------
#----------------------------------------
#Function to make most extreme values for the skewed tail positive, i.e., force all distributions to be right skewed, and order ICs by skewness.
rightskew=function(S,M=NULL,order.skew=TRUE) {
  #S: n x d matrix
  #A: d x d` corresponding to X = S A
  #If order = TRUE, then ICs are organized from HIGHEST to LOWEST skewness where skewness is forced to be positive for all ICs.
  nObs <- nrow(S)
  if(ncol(S)>nObs) stop('S must be n x d')
  skewness<-function(x,n = nObs) (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  skew=apply(S,2,skewness)
  sign=-1*(skew<0)+1*(skew>0)
  S.new=S%*%diag(sign)
  if(!is.null(M)) M.new=t(diag(sign))%*%M
  if(order.skew==TRUE) {
    skew.new=apply(S.new,2,skewness)
    perm=order(skew.new,decreasing=TRUE)
    S.new=S.new[,perm]
    if(!is.null(M)) M.new=M.new[perm,]
  }
  if(is.null(M)) {
    S.new
  } else
    return(list(S=S.new,M=M.new))
}
#---------------------------------------------
