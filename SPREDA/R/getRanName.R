getRanName <-
function(obj){
  if(obj[[1]]!="~") print("Please use the correct formula symbol")
  random=obj[[2]]
  m=length(random)
  if(m!=3) print("The length of the randome effects are not 3.")
  if(random[[1]]!="|") print("There is no conditional term in Random effects")
  cov.names=getnames(random[[2]])
  cond.names=getnames(random[[3]])
  res=list(lhs=cov.names, rhs=cond.names)
  return(res)
}
