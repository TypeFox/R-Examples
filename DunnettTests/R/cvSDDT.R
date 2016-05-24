cvSDDT <-
function(k,alpha=0.05,alternative="U",df=Inf,corr=0.5,corr.matrix=NA){
  if(k > 16){
    stop("The funtion is not applicable to the situations where the number of tests exceeds 16")
  }
  cvSet <- NULL
  ##calculate C1
  if(alternative=="U"){
    C1 <- qt(1-alpha,df)
    }
  if(alternative=="B"){
    C1 <- qt(1-alpha/2,df)
  }
  cvSet <- c(cvSet,C1)

  if(is.matrix(corr.matrix)==TRUE){
  Cj_SD.fun <- function(Cj_SD,j){
  	#find all the subsets of test statistics with subset size be i
     subset <- combn(1:k,j)
     ProbSet <- NULL
     if(alternative=="U"){
     	if(df==Inf){
         #multivariate normal distribution
          for(i in 1:ncol(subset)){
          	si <- subset[,i]
          	corr.matrix.si <- corr.matrix[si,si]
          	ProbSet[i]<- pmvnorm(lower=rep(-Inf,j),upper=rep(Cj_SD,j),corr=corr.matrix.si)
          	}
          }
        if(df!=Inf){
         #multivariate t distribution
          for(i in 1:ncol(subset)){
          	si <- subset[,i]
          	corr.matrix.si <- corr.matrix[si,si]
          	ProbSet[i]<- pmvt(lower=rep(-Inf,j),upper=rep(Cj_SD,j),df=df,corr=corr.matrix.si)
          	}
          }
        }
    if(alternative=="B"){
    if(df==Inf){
    	for(i in 1:ncol(subset)){
          	si <- subset[,i]
          	corr.matrix.si <- corr.matrix[si,si]
          	ProbSet[i]<- pmvnorm(lower=rep(-Cj_SD,j),upper=rep(Cj_SD,j),corr=corr.matrix.si)
          	} 
          }
    if(df!=Inf){
    	for(i in 1:ncol(subset)){
          	si <- subset[,i]
          	corr.matrix.si <- corr.matrix[si,si]
          	ProbSet[i]<- pmvt(lower=rep(-Cj_SD,j),upper=rep(Cj_SD,j),df=df,corr=corr.matrix.si)
          	}   
          }
        }
  min(ProbSet)-(1-alpha)
  }
  for(j in 2:k){
    Cj_SD <- uniroot(Cj_SD.fun,j=j,interval= c(0,10),tol=1e-05)$root
    cvSet <- c(cvSet,Cj_SD)
  }
  }

  if(is.matrix(corr.matrix)==FALSE){
  corr.matrix <- diag(1,k)
  corr.matrix[which(corr.matrix!=1)] <- corr
  Cj_SD.fun <- function(Cj_SD,j){
  if(alternative=="U"){
    if(df==Inf){
      #multivariate normal distribution
      Prob <- pmvnorm(lower=rep(-Inf,j),upper=rep(Cj_SD,j),corr=corr.matrix[1:j,1:j])-(1-alpha)  
    }
    if(df!=Inf){
      #multivariate t distribution
      Prob <- pmvt(lower=rep(-Inf,j),upper=rep(Cj_SD,j),df=df,corr=corr.matrix[1:j,1:j])-(1-alpha)  
    }
  }
  if(alternative=="B"){
    if(df==Inf){
      Prob <- pmvnorm(lower=rep(-Cj_SD,j),upper=rep(Cj_SD,j),corr=corr.matrix[1:j,1:j])-(1-alpha)
    }
    if(df!=Inf){
      Prob <- pmvt(lower=rep(-Cj_SD,j),upper=rep(Cj_SD,j),df=df,corr=corr.matrix[1:j,1:j])-(1-alpha)  
    }
  }
  Prob
  }
  for(j in 2:k){
    Cj_SD <- uniroot(Cj_SD.fun,j=j,interval=c(0,10),tol=1e-05)$root
    cvSet <- c(cvSet,Cj_SD)
  }
  }
  names(cvSet)=paste("c",seq(1,k,by=1),sep="")
  return(round(cvSet,3))
}
