###
### R routines for the R package mvmeta (c) Antonio Gasparrini 2014
#
checkPD <- 
function(mat, k, label="matrix") {
#
################################################################################
# 
  if(any(dim(mat)!=k)) stop(paste("wrong dimensions for '",label,"'",sep=""))
  eig <- eigen(mat)
  if(any(eig$values<0)) 
    stop(paste("values for '",label,
      "' do not produce a positive-definite matrix",sep=""))
#
  mat
}
