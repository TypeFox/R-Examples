coP_nrm <-
function(pitemNRM,Km)
{ # internal function to compute the P´s for NRM

  LAM <- matrix(pitemNRM,nrow=2,byrow=T)
  
  Z <- Km %*% LAM
  ez <- exp(Z)
  ezrs <- rowSums(ez)        
  ZQs <- ez / ezrs
  
  ZQs
}
