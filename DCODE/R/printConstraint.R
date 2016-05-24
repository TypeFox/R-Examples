printConstraint <-
function(l,nPep,pepList,nbAA=20){
  
  res=paste(paste(pepList[[nPep]][l[which(l<=nbAA^nPep)]],collapse=" + "), " = ",paste(tolower(pepList[[nPep]][l[which(l>nbAA^nPep)]-nbAA^nPep]), collapse=" + "))
  
  if(length(l)==1){
    if(l<=nbAA^nPep){
      res= paste(pepList[[nPep]][l[which(l<=nbAA^nPep)]],"= 0")
    }else{
      res=paste("0 =",tolower(pepList[[nPep]][l[which(l>nbAA^nPep)]-nbAA^nPep]))
    }
  }
  return(res)
}
