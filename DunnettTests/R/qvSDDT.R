qvSDDT <-
function(teststats,alternative="U",df=Inf,corr = 0.5,corr.matrix = NA){
  k <- length(teststats)
  if(k > 16){
    stop("The funtion is not applicable to situations where the number of tests exceeds 16")
  }
  if(is.matrix(corr.matrix)==FALSE){
  corr.matrix <- diag(1,k)
  corr.matrix[which(corr.matrix!=1)] <- corr
  }
  names(teststats) <- paste("H",1:k,sep="")
  
  #order the test statistcis
  if(alternative=="U"){od.index <- order(teststats,decreasing=TRUE)}
  if(alternative=="B"){od.index <- order(abs(teststats),decreasing=TRUE)}
  
  t.ordered <- teststats[od.index]
  corr.matrix <- corr.matrix[od.index,od.index]

  pvSet <- NULL #to store the unadjusted P values
  qvSet <- NULL #to store the adjusted P values
  
  for(j in 1:(k-1)){
    if(alternative=="U"){
    if(df==Inf){
      #multivariate normal distribution
      Pj <- 1-pmvnorm(lower=rep(-Inf,k-j+1),upper=rep(t.ordered[j],k-j+1),corr=corr.matrix[j:k,j:k])  
    }
    if(df!=Inf){
      #multivariate t distribution
      Pj <- 1-pmvt(lower=rep(-Inf,k-j+1),upper=rep(t.ordered[j],k-j+1),df=df,corr=corr.matrix[j:k,j:k]) 
    }
  }
    if(alternative=="B"){
    if(df==Inf){
      Pj <- 1-pmvnorm(lower=rep(-abs(t.ordered[j]),k-j+1),upper=rep(abs(t.ordered[j]),k-j+1),corr=corr.matrix[j:k,j:k])
    }
    if(df!=Inf){
      Pj <- 1-pmvt(lower=rep(-abs(t.ordered[j]),k-j+1),upper=rep(abs(t.ordered[j]),k-j+1),df=df,corr=corr.matrix[j:k,j:k])
    }
  }
  pvSet <- c(pvSet,Pj)
  qvSet <- c(qvSet,max(pvSet))
}

  ##calculate Pk
  if(alternative=="U"){
    Pk <- as.numeric(1-pt(t.ordered[k],df))
  }
  if(alternative=="B"){
    Pk <- as.numeric(2*(1-pt(abs(t.ordered[k]),df)))
  }
  pvSet <- c(pvSet,Pk)
  qvSet <- c(qvSet,max(pvSet))
  Results <- list(t.ordered,round(qvSet,3))
  names(Results) <- c("ordered test statistics","Adjusted P-values of ordered test statistics")
  return(Results)
}
