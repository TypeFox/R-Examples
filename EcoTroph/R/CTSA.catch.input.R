CTSA.catch.input=function(catch.group,smooth_type=NULL,sigmaLN_cst=NULL,pas=NULL,shift=NULL,smooth_param=NULL){
  nm=colnames(catch.group)
  fleet=nm[substring(nm,1,6)=='catch.'];n.fleet=length(fleet)
  catch=list()
  for(i in 1:n.fleet){catch[[fleet[i]]]=Transpose(create.smooth(catch.group,smooth_type=smooth_type,sigmaLN_cst=sigmaLN_cst,pas=pas,shift=shift,smooth_param=smooth_param),catch.group,fleet[i])}
  return(catch)
}