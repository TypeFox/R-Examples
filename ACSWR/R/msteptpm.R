msteptpm <-
function(TPM,m){
  if(m==1) return(TPM) else {
    temp <- TPM
    for(i in 1:(m-1)) temp=temp%*%TPM
    return(temp)
  }
}
