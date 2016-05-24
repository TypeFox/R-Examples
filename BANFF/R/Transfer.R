#####Transferring large statistical pvalue into testing statistics
Transfer=function(pvalue)
{
  rstat=-(qnorm(pvalue))
  return(rstat)
}
