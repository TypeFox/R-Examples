"recode" <-
function(x,oldcode=sort(unique(x)),newcode){
  if(length(oldcode)!=length(newcode))stop("The number of old and new codes do not match")
  newx<-x
  for(i in 1:length(oldcode)){
    newx[x==oldcode[i]]<-newcode[i]
  }
  return(newx)
}

