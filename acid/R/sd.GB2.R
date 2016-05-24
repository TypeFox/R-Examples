sd.GB2 <-
function(b,a,p,q){
  fmGB2<- km.GB2(b,a,p,q,k=1) # first moment
  smGB2<- km.GB2(b,a,p,q,k=2) # second moment
  varGB2  <- smGB2-fmGB2^2
  sdGB2   <- sqrt(varGB2)
  return(sdGB2)
}
