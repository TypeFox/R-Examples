# @title Variance-covariance matrix of a CUBE model with covariates
# @description Compute the variance-covariance matrix of parameter estimates of a CUBE model with covariates
#  for all the three parameters.
# @aliases varcovcubecov
# @usage varcovcubecov(m, ordinal, Y, W, Z, estbet, estgama, estalpha)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of covariates for explaining the uncertainty component
# @param W Matrix of covariates for explaining the feeling component
# @param Z Matrix of covariates for explaining the overdispersion component
# @param estbet Vector of the estimated parameters for the uncertainty component, with length equal to 
# NCOL(Y)+1 to account for an intercept term (first entry)
# @param estgama Vector of the estimated parameters for the  feeling component, with length equal to 
# NCOL(W)+1 to account for an intercept term (first entry)
# @param estalpha Vector of the estimated parameters for the overdispersion component, with length
#  equal to NCOL(Z)+1 to account for an intercept term (first entry)
# @details The function checks if the variance-covariance matrix is positive-definite: if not, 
# it returns a warning message and produces a matrix with NA entries.
# @seealso \code{\link{loglikCUBE}}, \code{\link{CUBE}} 
# @references
# Piccolo, D. (2014), Inferential issues on CUBE models with covariates, 
# \emph{Communications in Statistics - Theory and Methods}, \bold{44}, DOI: 10.1080/03610926.2013.821487 
#' @keywords internal

varcovcubecov <-
function(m,ordinal,Y,W,Z,estbet,estgama,estalpha){
  n<-length(ordinal)
  p<-NCOL(Y);q<-NCOL(W);v<-NCOL(Z);
  paivett<-logis(Y,estbet)
  csivett<-logis(W,estgama)
  phivett<-1/(-1+ 1/(logis(Z,estalpha)))
  # Probability
  probi<-paivett*(betabinomial(m,ordinal,csivett,phivett)-1/m)+1/m
  uui<-1-1/(m*probi)
  ubari<-uui+paivett*(1-uui)
  ### Matrix computations
  mats1<-auxmat(m,csivett,phivett,1,-1,1,0,1)
  mats2<-auxmat(m,csivett,phivett,0, 1,1,0,1)
  mats3<-auxmat(m,csivett,phivett,1,-1,1,1,1)
  mats4<-auxmat(m,csivett,phivett,0, 1,1,1,1)
  mats5<-auxmat(m,csivett,phivett,1, 0,1,1,1)
  ### for D_i vectors
  matd1<-auxmat(m,csivett,phivett,1,-1,2,0, 1)
  matd2<-auxmat(m,csivett,phivett,0, 1,2,0,-1)
  matd3<-auxmat(m,csivett,phivett,1,-1,2,1, 1)
  matd4<-auxmat(m,csivett,phivett,0, 1,2,1,-1)
  ### for H_i vectors
  math3<-auxmat(m,csivett,phivett,1,-1,2,2,-1)
  math4<-auxmat(m,csivett,phivett,0, 1,2,2,-1)
  math5<-auxmat(m,csivett,phivett,1, 0,2,2,-1)
  ###
  S1<-S2<-S3<-S4<-D1<-D2<-D3<-D4<-H3<-H4<-rep(NA,m);
  for(i in 1:n){
    S1[i]<-sum(mats1[1:ordinal[i],i])-mats1[ordinal[i],i]
    S2[i]<-sum(mats2[1:(m-ordinal[i]+1),i])-mats2[(m-ordinal[i]+1),i]
    S3[i]<-sum(mats3[1:ordinal[i],i])-mats3[ordinal[i],i]
    S4[i]<-sum(mats4[1:(m-ordinal[i]+1),i])-mats4[(m-ordinal[i]+1),i]
    D1[i]<-sum(matd1[1:ordinal[i],i])-matd1[ordinal[i],i]
    D2[i]<-sum(matd2[1:(m-ordinal[i]+1),i])-matd2[(m-ordinal[i]+1),i]
    D3[i]<-sum(matd3[1:ordinal[i],i])-matd3[ordinal[i],i]
    D4[i]<-sum(matd4[1:(m-ordinal[i]+1),i])-matd4[(m-ordinal[i]+1),i]
    H3[i]<-sum(math3[1:ordinal[i],i])-math3[ordinal[i],i]
    H4[i]<-sum(math4[1:(m-ordinal[i]+1),i])-math4[(m-ordinal[i]+1),i]
  }
  ###### Attention !!!
  S5<-colSums(mats5[1:(m-1),])
  H5<-colSums(math5[1:(m-1),])
  ##########
  CC<-S2-S1;   EE<-S3+S4-S5;
  DD<-D2-D1;   FF<-D3+D4;   GG<-H3+H4-H5;
  ### Computing vectors v and u
  vibe<-uui*(1-paivett)
  viga<-ubari*csivett*(1-csivett)*CC
  vial<-ubari*phivett*EE  #****
  ubebe<-uui*(1-paivett)*(1-2*paivett)
  ugabe<-ubari*csivett*(1-csivett)*(1-paivett)*CC
  ualbe<-ubari*phivett*(1-paivett)*EE  #****
  ugaga<-ubari*csivett*(1-csivett)*((1-2*csivett)*CC+csivett*(1-csivett)*(CC^2+DD))
  ualga<-ubari*phivett*csivett*(1-csivett)*(FF+CC*EE)   #****
  ualal<-ubari*phivett*(EE+phivett*(EE^2+GG))   #****
  ### Computing vectors g
  gbebe<-Hadprod(vibe,vibe)-ubebe                       
  ggabe<-Hadprod(viga,vibe)-ugabe                        
  galbe<-Hadprod(vial,vibe)-ualbe                        
  ggaga<-Hadprod(viga,viga)-ugaga                        
  galga<-Hadprod(vial,viga)-ualga
  galal<-Hadprod(vial,vial)-ualal
  ### Expanding matrices Y, W, Z
  YY<-cbind(1,Y); WW<-cbind(1,W); ZZ<-cbind(1,Z);
  ### Elements of the Information matrix
  infbebe<-t(YY)%*%(Hadprod(YY,gbebe))
  infgabe<-t(WW)%*%(Hadprod(YY,ggabe)) 
  infalbe<-t(ZZ)%*%(Hadprod(YY,galbe)) 
  infgaga<-t(WW)%*%(Hadprod(WW,ggaga))         
  infalga<-t(ZZ)%*%(Hadprod(WW,galga))                 
  infalal<-t(ZZ)%*%(Hadprod(ZZ,galal))            
  infbega<-t(infgabe)
  infbeal<-t(infalbe)
  infgaal<-t(infalga)
  ### Assembling Information matrix...
  #   # matinf<-rbind(
  #     cbind(infbebe,infbega,infbeal),
  #     cbind(infgabe,infgaga,infgaal),
  #     cbind(infalbe,infalga,infalal));
  ### Variance-covariance matrix
  npai<-NCOL(Y)+1
  ncsi<-NCOL(W)+1
  nphi<-NCOL(Z)+1
  nparam<-npai+ncsi+nphi
  matinf<-matrix(NA,nrow=nparam,ncol=nparam)
  for (i in 1:npai){
    matinf[i,] <- matrix(c(infbebe[i,],infbega[i,],infbeal[i,]),nrow=1,byrow=TRUE)
  }
  for (i in (npai+1):(npai+ncsi)){
    matinf[i,] <- matrix(c(infgabe[i-npai,],infgaga[i-npai,],infgaal[i-npai,]),nrow=1,byrow=TRUE)
  }
  for (i in (npai+ncsi+1):nparam){
    matinf[i,] <- matrix(c(infalbe[i-npai-ncsi,],infalga[i-npai-ncsi,],infalal[i-npai-ncsi,]),nrow=1,byrow=TRUE)
  }
  
  
  if(any(is.na(matinf))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=nparam,ncol=nparam)
  } else {
    if(det(matinf)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=nparam,ncol=nparam)
    } else {
      varmat<-solve(matinf)
    }
  }
  
  return(varmat)
}
