marg_param <- function(lev,type){

#        [C,M,G] = marg_param(lev,type)
#
# Creates matrices C and M for the marginal parametrization
# of the probability vector for a vector of categorical variables
# 
# INPUT:
# lev:  vector containing the number of levels of each variable
# type: vector with elements 'l', 'g', 'c', 'r' indicating the type of logit
#
# OUTPUT:
# C:    matrix of constrats (the first sum(lev)-length(r) elements are
#       referred to univariate logits)
# M:    marginalization matrix with elements 0 and 1
# G:    corresponding design matrix for the corresponding log-linear model

# preliminaries
  r = length(lev)
# create sets of parameters
  S = t(apply(sq(r,1),1,rev))
  if(r>=2) S = rbind(S,t(apply(sq(r,2),1,rev)))
# create matrices
  C = NULL; M = NULL; G = NULL
  for(i in 1:nrow(S)){
    si = S[i,]; 
    Ci = 1;  # to update matrix C
    for(h in 1:r){
      if(si[h]==1){
        I = diag(lev[h]-1)
        Ci = Ci%x%cbind(-I,I)
      }
    }
    if(i==1) C = Ci else C = blkdiag(C,Ci)
    Mi = 1;  # to update matrix M
    for(h in 1:r){
      if(si[h]==1){
        I = diag(lev[h]-1); T1 = matrix(1,lev[h]-1,lev[h]-1); T1[upper.tri(T1)]=0; ze = matrix(0,lev[h]-1,1);
        if(type[h]=="l") Mi = Mi%x%rbind(cbind(I,ze),cbind(ze,I))
        if(type[h]=="g") Mi = Mi%x%rbind(cbind(T1,ze),cbind(ze,t(T1)))
        if(type[h]=="c") Mi = Mi%x%rbind(cbind(I,ze),cbind(ze,t(T1)))
        if(type[h]=="r") Mi = Mi%x%rbind(cbind(T1,ze),cbind(ze,I))
      }else{
        Mi = Mi%x%matrix(1,1,lev[h])
      }
    }
    M = rbind(M,Mi)
    Gi = 1;  # for the desgin matrix
    for(h in 1:r){
      if(si[h]==1){
        Gi = Gi%x%t(cbind(-matrix(1,lev[h]-1,1),diag(lev[h]-1)))/2
      }else{
        Gi = Gi%x%matrix(1,lev[h],1)
      }
    }
    G = cbind(G,Gi)
  }
  out = list(C=C,M=M,G=G)
  out
}