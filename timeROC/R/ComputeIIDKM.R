# function to compute iid of KM estimator of the censoring survival distribution
# {{{ Input
#times : observed time
#status : 1 if non censored, 0 if censored
# }}}
# {{{ Output
#iid.mat : matrix of the iid representation of KM for all time in times
# }}}
Compute.iid.KM <- function(times,status){
  #browser()
  times <- times[order(times)]
  status <- status[order(times)] 
  n <- length(times)
  mat.data<-cbind(times,as.numeric(status==0))
  colnames(mat.data)<-c("T","indic.Cens")
  # compute the empirical survival function corresponding to the counting process 1(\tilde{eta}=0, \tilde{T}<=t)
  hatSdeltaCensTc<-1-cumsum(mat.data[,c("indic.Cens")])/n  
  # Build the matrix required for computing  dM_C(u) for all time u (all observed times \tilde{T}_i)
  temp1 <- cbind(mat.data[,c("T","indic.Cens")],1-(1:n)/n,hatSdeltaCensTc)
  temp1 <- rbind(c(0,0,1,1),temp1) # Add the first row corresponding to time t=0
  colnames(temp1)<-c("T","indic.Cens","hatSTc","hatSdeltaCensTc")
  # compute hazard function of the censoring
  lambdaC<-(temp1[-1,"indic.Cens"])/(n:1)  
  # Add the column of the hazard function of the censoring (equal to 0 at time t=0)
  temp1<-cbind(temp1,c(0,lambdaC))
  colnames(temp1)[ncol(temp1)]<-"lambdaC"
  # Cumulative hazard of censoring
  LambdaC<-cumsum(lambdaC)         
  # Add the column of the cumulative hazard function of the censoring (equal to 0 at time t=0)
  temp1 <- cbind(temp1,c(0,LambdaC))
  colnames(temp1)[ncol(temp1)]<-"LambdaC"
  temp2<-temp1[-1,]
  # compute  martingale of censoring \hat{M}_{C_i}(u) for all time u (all observed times \tilde{T}_i) using previous matrix
  # We obtain a matrix. Each column contains the vector of M_{C_i}(\tilde{T}_j) for  all j.
  hatMC<-matrix(NA,n,n)
  for (i in 1:n){
    hatMC[,i] <-temp2[i,2]*as.numeric(temp2[i,1]<=temp2[,"T"])- c(temp2[0:i,"LambdaC"], rep(temp2[i,6],(n-i)))
  }  
  # In order to draw martingale paths
  #matplot(mat.data[,"T"],hatMC,type="l")
  #lines(mat.data[,"T"],rowMeans(hatMC),lwd=5)  
  # Compute d \hat{M}_{C_i} (u) for all time u (all observed times \tilde{T}_i)
  dhatMC<-rbind(hatMC[1,],hatMC[-1,]-hatMC[-nrow(hatMC),])
  # Compute d \hat{M}_{C_i} (u)/(S_{\tilde{T}}(u)) for all time u (all observed times \tilde{T}_i)
  # We need this for integrals in the martingale representation of the Kaplan-Meier estimator of the censoring survival function
  # function to divide d \hat{M}_{C_i} (u) by (S_{\tilde{T}}(u))
  MulhatSTc<-function(v){
    n <- length(v)
    v/c(1,1-(1:(n-1))/n)      # c(1,1-(1:(n-1))/n) is the at risk probability (S_{\tilde{T}}(u))
  }
  # apply the function for each column (corresponding to the
  # vector M_{C_i}(u)  for all time u (all observed times \tilde{T}_i), 
  # time \tilde{T}_i corresponds to the i-th row of the matrix)
  dhatMCdivST<-apply(dhatMC,2,MulhatSTc)
  # Compute \int_0^{\tilde{T}_j} d{ \hat{M}_{C_l} (u) } / (S_{\tilde{T}}(u)) for each subject l, we compute for all time \tilde{T}_j.
  # l=column, j=row
  MatInt0TcidhatMCksurEff<-apply(dhatMCdivST,2,cumsum)  # (Remark : on of the row corresponds to the previous step...) 
  colnames(MatInt0TcidhatMCksurEff)<-paste("M_{C_",1:length(times),"}",sep="")
  rownames(MatInt0TcidhatMCksurEff)<-times  
  return(MatInt0TcidhatMCksurEff)  
}

