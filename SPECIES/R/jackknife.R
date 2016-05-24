####################################################################################################

##################################################################################################
jackknife= function(n,k=5,conf=0.95){
  
  ## ==============================================================================================
  ## Purpose: This function calculates the Jackknife estimator  by Burnham and Overton 1978 and 1979.
  ## input:   n --- frequency and frequency data n, matrix format, numeral data frame
  ##          k --- integer, jackknife order
  ##          conf--- confidence level. The order was automatically selected using the test at the
  ##                  significance level=1-conf. In the original paper, the authors have used significance
  ##                  level=0.05. The confidence  interval at the selected or specified order is calculated
  ##                  at this confidence level.
  ## output:  jackknife estimator, order, and standard error and confidence interval.
  ## ===============================================================================================


  if (k!=round(k)||k<0) stop("Error: The Jackknife order k must be an positive integer!")
  if(is.numeric(conf)==FALSE||conf>1||conf<0) stop("Error: confidence level must be a numerical value between 0 and 1, e.g. 0.95")
  if (ncol(n)!=2||is.numeric(n[,1])==FALSE||is.numeric(n[,2])==FALSE) {
    stop("Error: The frequency of frequencies data n must be a matrix.")
  }

  k0=k
  k=min(nrow(n)-1,10)
  n=as.matrix(n)
  m=max(n[,1])
  ntemp=cbind(c(1:m),rep(0,m))
  ntemp[n[,1],2]=n[,2]
  n=ntemp
  colnames(n)=NULL
  rownames(n)=NULL 
  cat("\n") 
  
  if(k0>k){
    cat("Warning: the Jackknife order in your data can not be larger than",k,"!","\n")
    cat("Up to orde", k, "jackknife estimate will be estiamted.","\n\n")
  }
  total = sum(n[,2])
  gene = matrix(0,k+1,5)
  gene[1,1] = total
  gene[1,2] = 0
  gene[1,5] = 0
  

  ## point estimator (First Eq. on page 935, Burnham and Overton 1979.
  ## Testing procedure( last and last second Equations on page 929, Burnham and Overton 1979)
  ## Standard error (First Eq. on page 930, Burnham and Overton 1979)
  
  for (i in 1:k){
    gene[i+1,1] = total
    gene[i+1,4] = total
    for (j in 1:i){
      gene[i+1,1] = gene[i+1,1]+(-1)^(j+1)*2^i*dbinom(j,i,.5)*n[j,2]
      gene[i+1,4] = gene[i+1,4]+(-1)^(j+1)*2^i*dbinom(j,i,.5)*n[j,2]*prod(1:j)
    }
    gene[i+1,2] = -gene[i+1,1]
    for (j in 1:i){
      gene[i+1,2] =gene[i+1,2] +((-1)^(j+1)*2^i*dbinom(j,i,.5)+1)^2*n[j,2]
    }
    gene[i+1,2]=gene[i+1,2]+sum(n[(i+1):nrow(n),2])
    gene[i+1,2] = sqrt(gene[i+1,2])
  }
  if(k>1){
    for (i in 2:(k)){
      gene[i,3] = -(gene[i+1,1] -gene[i,1])^2/(total-1)
      for (j in 1:(i-1)){
        gene[i,3] =gene[i,3] +((-1)^(j+1)*2^(i)*dbinom(j,i,.5)-(-1)^(j+1)*2^(i-1)*dbinom(j,i-1,.5))^2*n[j,2]*total/(total-1)
      }
      gene[i,3] = gene[i,3]+n[i,2]*total/(total-1)
      #gene[i,3] = gene[i,3]+((-1)^(i)*2^(i)*dbinom(i,i,.5))^2*n[i,2]*total/(total-1)
      gene[i,3] = sqrt(gene[i,3])
      gene[i,5] = (gene[i+1,1]-gene[i,1])/gene[i,3]	    
    }
  }
  coe=qnorm(1-(1-conf)/2,0,1)
  x=(gene[2:(k+1),5]<coe)
  if (length(x[x==TRUE])==0) {
    jackest = gene[k+1,1]
    sej= gene[k+1,2]
  }  else{
    indicator=c(2:(k+1))
    if(length(x[x==TRUE])>0) {
      jackest = gene[indicator[x][1],1]
      sej= gene[indicator[x][1],2]
      order=c(1:(k+1))[indicator[x][1]]-1
    }
  }
  
  if(k0<= order){
    jackest=gene[k0+1,1]
    sej= gene[k0+1,2]
    order=k0
  }else if(order<k0){
    cat("Your specified order is larger than that determined by the test,","\n")
    cat("Therefore the order from the test is used.","\n\n")
  }
  CI0=matrix(c(jackest-coe*sej,jackest+coe*sej),1,2)
  colnames(CI0)=c("lb","ub")
  return(list(JackknifeOrder=order, Nhat=round(jackest),SE=sej,CI=round(CI0)))
}

