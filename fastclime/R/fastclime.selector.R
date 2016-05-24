# fastclime.lambda(): Function used to select the solution path for a given lambda         #
# Authors: Haotian Pang, Han Liu and Robert Vanderbei                                      #
# Emails: <hpang@princeton.edu>, <hanliu@princeton.edu> and <rvdb@princetonedu>            #
# Date: April 22th 2016                                                           #
# Version: 1.4.1					                                                       #
#------------------------------------------------------------------------------------------#

fastclime.selector <- function(lambdamtx, icovlist, lambda)
{

  gcinfo(FALSE)
  d<-dim(icovlist[[1]])[2]
  maxnlambda<-dim(lambdamtx)[1]
  icov<-matrix(0,d,d)
  adaj<-matrix(0,d,d)
  seq<-rep(0,d)
  threshold<-1e-5
  status<-0
  
  for(i in 1:d)
  {
    
      temp_lambda<-which(lambdamtx[,i]>lambda)
      seq[i]<-length(temp_lambda)
      
      if((seq[i]+1)>maxnlambda)
      {
        status<-1
        icov[,i]<-icovlist[[seq[i]]][,i]
      }
      else{
        icov[,i]<-icovlist[[seq[i]+1]][,i]
      }
     
  }
  
  icov <- icov*(abs(icov)<= abs(t(icov)))+ t(icov)*(abs(icov)> abs(t(icov)))

  tmpicov<-icov
  diag(tmpicov)<-0
  adaj<-Matrix(tmpicov>threshold, sparse=TRUE)*1
 

  
  sparsity<-(sum(adaj))/(d^2-d)
  
  if(status==1)
  {
    cat("Some columns do not reach the required lambda!\n You may want to increase lambda.min or use a larger nlambda. \n")
  }
  
  rm(temp_lambda,seq,d,threshold)
  gc()
  
  result<-list("icov"=icov, "adaj"=adaj,"sparsity"=sparsity)
  class(result)="fastclime.selector"
  
  return(result)


}


