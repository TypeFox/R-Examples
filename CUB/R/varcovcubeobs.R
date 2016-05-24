# @title Variance-covariance matrix for CUBE models based on the observed information matrix
# @description Compute the variance-covariance matrix of parameter estimates for a CUBE model without covariates 
# as the inverse of the observed information matrix.
# @aliases varcovcubeobs
# @usage varcovcubeobs(m, pai, csi, phi, freq)
# @param m Number of ordinal categories
# @param pai Uncertainty parameter
# @param csi Feeling parameter
# @param phi Overdispersion parameter
# @param freq Vector of the observed absolute frequencies
# @details The function checks if the variance-covariance matrix is positive-definite: if not, 
# it returns a warning message and produces a matrix with NA entries.
# @seealso \code{\link{cube000}}, \code{\link{varcovcubeexp}}
# @references
# Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
# \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786
#' @keywords internal


varcovcubeobs <-
function(m,pai,csi,phi,freq){
  pr<-probcube(m,pai,csi,phi)
  sum1<-sum2<-sum3<-sum4<-rep(NA,m);### sums computation
  d1<-d2<-h1<-h2<-h3<-h4<-rep(NA,m); ### derivatives computation
  ### sum1; sum2; sum3; sum4; sum5; as in Iannario (2013), "Comm. in Stat. Theory & Methods"
  
  for(jr in 1:m){
    seq1<-1/((1-csi)+phi*((1:jr)-1))                ### a(k)
    seq2<-1/((csi)+phi*((1:(m-jr+1))-1))            ### b(k)
    seq3<-((1:jr)-1)/((1-csi)+phi*((1:jr)-1))
    seq4<-((1:(m-jr+1))-1)/((csi)+phi*((1:(m-jr+1))-1))
    dseq1<-seq1^2
    dseq2<-seq2^2 
    hseq1<-dseq1*((1:jr)-1)
    hseq2<-dseq2*((1:(m-jr+1))-1)
    hseq3<-dseq1*((1:jr)-1)^2
    hseq4<-dseq2*((1:(m-jr+1))-1)^2
    #############
    sum1[jr]<-sum(seq1)-seq1[jr]
    sum2[jr]<-sum(seq2)-seq2[m-jr+1]
    sum3[jr]<-sum(seq3)-seq3[jr]
    sum4[jr]<-sum(seq4)-seq4[m-jr+1]
    d1[jr]<- sum(dseq1)-dseq1[jr]
    d2[jr]<- -(sum(dseq2)-dseq2[m-jr+1])
    h1[jr]<- -(sum(hseq1)-hseq1[jr])
    h2[jr]<- -(sum(hseq2)-hseq2[m-jr+1])
    h3[jr]<- -(sum(hseq3)-hseq3[jr])
    h4[jr]<- -(sum(hseq4)-hseq4[m-jr+1])
  }
  seq5<-(0:(m-2))/(1+phi*(0:(m-2)))           ### c(k) correction ?!
  sum5<-sum(seq5)                             ### correction ?!
  h5<- -sum(seq5^2)                            ### correction ?!
  
  ### Symbols as in Iannario (2013), "Comm. in Stat.", ibidem (DP notes)
  uuur<-1-1/(m*pr)
  ubar<-uuur+pai*(1-uuur)
  vbar<-ubar-1
  aaar<-sum2-sum1
  bbbr<-sum3+sum4-sum5
  cccr<-h3+h4-h5
  dddr<-h2-h1
  eeer<-d2-d1
  ###### dummy product
  prodo<-freq*ubar
  ######
  infpaipai<-sum(freq*uuur^2)/pai^2
  infpaicsi<-sum(prodo*(uuur-1)*aaar)/pai
  infpaiphi<-sum(prodo*(uuur-1)*bbbr)/pai
  infcsicsi<-sum(prodo*(vbar*aaar^2-eeer))
  infcsiphi<-sum(prodo*(vbar*aaar*bbbr-dddr))
  infphiphi<-sum(prodo*(vbar*bbbr^2-cccr))
  ### Information matrix 
  inform<-matrix(NA,nrow=3,ncol=3)
  inform[1,1]<-infpaipai;   inform[1,2]<-infpaicsi;   inform[1,3]<-infpaiphi;
  inform[2,1]<-infpaicsi;   inform[2,2]<-infcsicsi;   inform[2,3]<-infcsiphi;
  inform[3,1]<-infpaiphi;   inform[3,2]<-infcsiphi;   inform[3,3]<-infphiphi;
  ### Var-covar matrix 
  nparam<-3
  if(any(is.na(inform))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=nparam,ncol=nparam)
  } else {
    if(det(inform)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=nparam,ncol=nparam)
    } else {
      varmat<-solve(inform)
    }
  }
  
  return(varmat)
}
